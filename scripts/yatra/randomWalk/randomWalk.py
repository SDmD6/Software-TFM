#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import networkx as nx
from scipy import sparse

__all__ = ["load_network_tsv", "build_graph_and_WT", "rwr", "run_rwwr"]

def load_network_tsv(path):
    """
    Lee una red en TSV aceptando dos variantes:
    1) Con cabecera: (from,to,score)  o  (u,v,w)
    2) Sin cabecera: se interpretan como (u,v,w)

    Devuelve un DataFrame con columnas ["u","v","w"], simétrico, sin bucles y con pesos max por par.
    """
    try:
        df = pd.read_csv(path, sep="\t")
        if {"from","to","score"}.issubset(df.columns):
            df = df.rename(columns={"from":"u","to":"v","score":"w"})
        elif {"u","v","w"}.issubset(df.columns):
            pass
        else:
            # Forzamos lectura sin cabecera
            df = pd.read_csv(path, sep="\t", header=None, names=["u","v","w"])
    except Exception:
        # Fallback sin cabecera
        df = pd.read_csv(path, sep="\t", header=None, names=["u","v","w"])

    # Limpieza y tipado
    df = df.dropna()
    df["w"] = pd.to_numeric(df["w"], errors="coerce").fillna(0.0)

    # Simetría (red no dirigida) y depuración
    df_u = df[["u","v","w"]]
    df_v = df.rename(columns={"u":"v","v":"u"})[["u","v","w"]]
    df   = pd.concat([df_u, df_v], ignore_index=True)
    df   = df[df["u"] != df["v"]]  # sin bucles
    df   = df.groupby(["u","v"], as_index=False)["w"].max()  # max peso por par
    return df

def build_graph_and_WT(df_edges):
    """
    Construye un grafo no dirigido y devuelve:
    - nodes: lista de nodos en orden
    - WT: matriz dispersa (W^T) con W fila-estocástica (normalizada por filas)
    """
    G = nx.from_pandas_edgelist(df_edges, "u", "v", edge_attr="w", create_using=nx.Graph())
    nodes = list(G.nodes())
    if len(nodes) == 0:
        raise RuntimeError("La red está vacía tras la limpieza/simetrización.")

    A = nx.to_scipy_sparse_array(G, nodelist=nodes, weight="w", dtype=float, format="csr")
    # Normalización por filas
    degrees = np.asarray(A.sum(axis=1)).ravel()
    degrees[degrees == 0.0] = 1.0
    D_inv = sparse.diags(1.0 / degrees)
    W = D_inv @ A
    WT = W.transpose().tocsr()
    return nodes, WT

def rwr(WT, p0, alpha=0.75, tol=1e-6, max_iter=10000):
    """
    Random Walk with Restart con matriz fila-estocástica (usando W^T).
    p_{t+1} = alpha * (W^T * p_t) + (1 - alpha) * p0
    """
    p = p0.copy()
    for _ in range(max_iter):
        p_new = alpha * (WT @ p) + (1.0 - alpha) * p0
        if np.linalg.norm(p_new - p, 1) < tol:
            return p_new
        p = p_new
    return p  # si no converge, devuelve el último (poco habitual con estos parámetros)

def run_rwwr(network_file, training_file, output_file,
             alpha=0.75, tol=1e-6, max_iter=10000, drop_seed_scores=False):
    """
    Ejecuta el flujo clásico:
    - Lee y limpia la red.
    - Construye WT.
    - Construye p0 uniforme en semillas presentes.
    - Ejecuta RWR.
    - Escribe TSV (gene, score) ordenado.
    """
    # Red
    net = load_network_tsv(network_file)
    nodes, WT = build_graph_and_WT(net)

    # Semillas
    with open(training_file, "r", encoding="utf-8") as f:
        seeds_all = [l.strip() for l in f if l.strip()]
    seeds = [g for g in seeds_all if g in nodes]
    if len(seeds) == 0:
        raise RuntimeError("Ningún gen de entrenamiento está en la red.")

    idx_train = [nodes.index(g) for g in seeds]
    p0 = np.zeros(len(nodes), dtype=float)
    p0[idx_train] = 1.0 / len(idx_train)

    p = rwr(WT, p0, alpha=alpha, tol=tol, max_iter=max_iter)
    scores = np.asarray(p).ravel()

    if drop_seed_scores:
        mask = np.zeros(len(nodes), dtype=float)
        mask[idx_train] = 1.0
        scores = scores * (1.0 - mask)

    out = pd.DataFrame({"gene": nodes, "score": scores}).sort_values("score", ascending=False)
    out.to_csv(output_file, sep="\t", index=False)
