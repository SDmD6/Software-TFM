# scripts/run_random_walk.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
})

# Ruta a tu ejecutable de Python (ajusta si usas venv)
python_path <- "C:/Users/Usuario/AppData/Local/Programs/Python/Python313/python.exe"
yatra_py     <- "scripts/yatra/randomWalk/run_random_walk.py"

# Lee par de enfermedades
disease_table <- read_tsv("data/disease_pairs.tsv", show_col_types = FALSE)

# Parámetros del RWR
alpha   <- 0.75
tol     <- 1e-6
maxiter <- 10000
drop_seeds <- FALSE  # pon TRUE si quieres score=0 en semillas

for (disease in disease_table$Disease) {
  
  # Carpetas
  prio_dir <- glue("output/{disease}/prioritization")
  if (!dir.exists(prio_dir)) dir.create(prio_dir, recursive = TRUE)
  
  # Entradas
  # Si guardas la red física por enfermedad:
  net_file <- glue("output/network/{disease}_network_physical.tsv")
  # Si por el contrario la tenías global en output/network/, usa:
  # net_file <- glue("output/network/{disease}_network_physical.tsv")
  
  # Normaliza nombres de columnas -> TSV para Yatra
  net_yatra <- glue("{prio_dir}/network.tsv")
  read_tsv(net_file, show_col_types = FALSE) %>%
    rename(u = from_gene, v = to_gene, w = combined_score) %>%
    write_tsv(net_yatra)
  
  # Semillas
  gene_file  <- glue("data/training_genes/{disease}_training_genes.tsv")
  training   <- glue("{prio_dir}/training.txt")
  read_tsv(gene_file, show_col_types = FALSE) %>%
    pull(gene) %>%
    write_lines(training)
  
  # Salidas
  out_classic <- glue("{prio_dir}/ranking_classic.tsv")
  out_expr    <- glue("{prio_dir}/ranking_expression.tsv")
  expr_file   <- glue("{prio_dir}/expression.tsv")  # opcional
  
  # Ejecuta CLÁSICO
  cmd_classic <- glue('"{python_path}" "{yatra_py}" ',
                      '--network "{net_yatra}" --train "{training}" --out "{out_classic}" ',
                      "--mode classic ",
                      "--alpha {alpha} --tol {tol} --max-iter {maxiter} ",
                      if (drop_seeds) "--drop-seed-scores" else "")
  message("[Yatra-classic] ", disease)
  system(cmd_classic, ignore.stdout = FALSE, ignore.stderr = FALSE)
  
  # Ejecuta EXPRESSION (si hay expression.tsv)
  if (file.exists(expr_file)) {
    cmd_expr <- glue('"{python_path}" "{yatra_py}" ',
                     '--network "{net_yatra}" --train "{training}" --expression "{expr_file}" --out "{out_expr}" ',
                     "--mode expression ",
                     "--alpha {alpha} --tol {tol} --max-iter {maxiter} ",
                     if (drop_seeds) "--drop-seed-scores" else "")
    message("[Yatra-expression] ", disease)
    system(cmd_expr, ignore.stdout = FALSE, ignore.stderr = FALSE)
  }
}