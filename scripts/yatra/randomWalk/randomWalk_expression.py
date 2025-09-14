#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from randomWalk import load_network_tsv, build_graph_and_WT, rwr

__all__ = ["run_rwwr_expression"]

def _read_expression(path):
    df = pd.read_csv(path, sep="\t")
    if not {"gene","value"}.issubset(df.columns):
        raise RuntimeError("expression.tsv debe tener columnas: gene, value")
    df["value"] = pd.to_numeric(df["value"], errors="coerce").fillna(0.0)
    df.loc[df["value"] < 0, "value"] = 0.0
    return dict(zip(df["gene"], df["value"]))

def run_rwwr_expression(network_file, training_file, expression_file, output_file,
                        alpha=0.75, tol=1e-6, max_iter=10000, drop_seed_scores=False):
    """
    Igual que run_rwwr, pero p0 se pondera por la expresión de las semillas.
    """
    # Red
    net = load_network_tsv(network_file)
    nodes, WT = build_graph_and_WT(net)

    # Semillas
    with open(training_file, "r", encoding="utf-8") as f:
        seeds_all = [l.strip() for l in f if l.strip()]

    expr = _read_expression(expression_file)

    present, values = [], []
    for g in seeds_all:
        if g in nodes:
            present.append(g)
            values.append(float(expr.get(g, 0.0)))

    if len(present) == 0:
        raise RuntimeError("Ningún gen de entrenamiento está en la red.")
    total = sum(values)
    if total <= 0:
        raise RuntimeError("Todas las semillas tienen expresión 0 o NA en expression.tsv.")

    # p0 ponderado por expresión (normalizado a 1)
    p0 = np.zeros(len(nodes), dtype=float)
    for g, v in zip(present, values):
        p0[nodes.index(g)] = v / total

    p = rwr(WT, p0, alpha=alpha, tol=tol, max_iter=max_iter)
    scores = np.asarray(p).ravel()

    if drop_seed_scores:
        mask = np.zeros(len(nodes), dtype=float)
        for g in present:
            mask[nodes.index(g)] = 1.0
        scores = scores * (1.0 - mask)

    out = pd.DataFrame({"gene": nodes, "score": scores}).sort_values("score", ascending=False)
    out.to_csv(output_file, sep="\t", index=False)
