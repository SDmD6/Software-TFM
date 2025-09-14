### =================== PRIORIZACIÓN CON YATRA (clásico o expresión) =================== ###
suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
})

# --- Configuración ---
python_path  <- "C:/Users/Usuario/AppData/Local/Programs/Python/Python313/python.exe"
yatra_script <- "scripts/yatra/randomWalk/run_random_walk.py"

# Parámetros RWR (se pasarán si tu script Python los soporta)
alpha     <- 0.75
tol       <- 1e-6
max_iter  <- 10000
drop_seeds <- FALSE  # TRUE si quieres score=0 para semillas (si tu Python lo soporta)

# --- Utilidades internas ---
first_existing <- function(paths) {
  hits <- paths[file.exists(paths)]
  if (length(hits) > 0) hits[1] else character(0)
}

# --- Enfermedades ---
disease_table <- read_tsv("data/disease_pairs.tsv", show_col_types = FALSE)
stopifnot("Disease" %in% names(disease_table))

for (disease in disease_table$Disease) {
  
  # === Entradas esperadas ===
  # Red física: probamos varias convenciones de ruta
  net_file <- first_existing(c(
    glue("output/{disease}/network/{disease}_network_physical.tsv"),
    glue("output/{disease}/network/network_physical.tsv"),
    glue("output/network/{disease}_network_physical.tsv"),
    glue("output/network/{disease}_physical.tsv"),
    glue("output/network/{disease}.tsv")
  ))
  if (length(net_file) == 0) {
    stop(glue("No existe la red física para {disease}. Ejecuta scripts/build_network.R.\nProbadas:\n- {paste(c(
      glue('output/{disease}/network/{disease}_network_physical.tsv'),
      glue('output/{disease}/network/network_physical.tsv'),
      glue('output/network/{disease}_network_physical.tsv'),
      glue('output/network/{disease}_physical.tsv'),
      glue('output/network/{disease}.tsv')
    ), collapse = '\n- ')}"))
  }
  
  train_file <- glue("data/training_genes/{disease}_training_genes.tsv")
  if (!file.exists(train_file)) {
    stop(glue("No existe el training para {disease}: {train_file}\nEjecuta scripts/read_input.R antes."))
  }
  
  expression_file <- glue("data/expression/{disease}_expression.tsv")  # opcional
  
  # === Carpeta de salida ===
  yatra_dir <- glue("output/{disease}/prioritization")
  if (!dir.exists(yatra_dir)) dir.create(yatra_dir, recursive = TRUE)
  
  # === Preparar red para Yatra (from, to, score) ===
  net_df <- read_tsv(net_file, show_col_types = FALSE)
  
  cand_from  <- intersect(c("from_gene","from","source","Source","u","protein1","node1"), names(net_df))
  cand_to    <- intersect(c("to_gene","to","target","Target","v","protein2","node2"),   names(net_df))
  cand_score <- intersect(c("combined_score","score","weight","Weight","w","combined_score_scaled"), names(net_df))
  
  if (length(cand_from) == 0 || length(cand_to) == 0 || length(cand_score) == 0) {
    stop("No encuentro columnas equivalentes a (from, to, score) en: ", net_file,
         "\nColumnas disponibles: ", paste(names(net_df), collapse = ", "))
  }
  
  net_df <- net_df %>%
    rename(from = !!sym(cand_from[1]),
           to   = !!sym(cand_to[1]),
           score= !!sym(cand_score[1])) %>%
    select(from, to, score) %>%
    mutate(from = as.character(from),
           to   = as.character(to),
           score= suppressWarnings(as.numeric(score))) %>%
    filter(!is.na(from), !is.na(to), !is.na(score), from != to)
  
  # eliminar duplicados si existieran (conservar el máximo score)
  net_df <- net_df %>%
    group_by(from, to) %>%
    summarise(score = max(score), .groups = "drop")
  
  net_yatra <- glue("{yatra_dir}/network.tsv")
  write_tsv(net_df, net_yatra)
  
  # === Preparar lista de genes (uno por línea) ===
  genes_tbl <- read_tsv(train_file, show_col_types = FALSE)
  if (!"gene" %in% names(genes_tbl)) names(genes_tbl)[1] <- "gene"
  
  genes <- genes_tbl[["gene"]]
  genes <- as.character(genes)
  genes <- unique(genes[!is.na(genes) & genes != ""])
  
  if (length(genes) == 0) {
    warning(glue("Training vacío tras limpieza para {disease}. Se omite."))
    next
  }
  
  genes_yatra <- glue("{yatra_dir}/training.txt")
  write_lines(genes, genes_yatra)
  
  # === Detectar modo ===
  modo <- if (file.exists(expression_file)) "expression" else "classic"
  yatra_out <- glue("{yatra_dir}/yatra_ranking.tsv")
  
  # === Comando Python ===
  base_cmd <- c(
    shQuote(yatra_script),
    "--network", shQuote(net_yatra),
    "--train",   shQuote(genes_yatra),
    "--out",     shQuote(yatra_out)
  )
  # Solo pasa hiperparámetros si tu run_random_walk.py los soporta
  base_cmd <- c(base_cmd,
                "--alpha",   as.character(alpha),
                "--tol",     as.character(tol),
                "--max-iter",as.character(max_iter))
  if (drop_seeds) base_cmd <- c(base_cmd, "--drop-seed-scores")
  
  if (modo == "expression") {
    base_cmd <- c(base_cmd, "--mode", "expression", "--expression", shQuote(expression_file))
  } else {
    base_cmd <- c(base_cmd, "--mode", "classic")
  }
  
  message("[Yatra-", modo, "] ", disease)
  message("[INFO] Ejecutando modo ", modo, "...")
  # Capturamos stdout/err para debug
  py_out <- tryCatch(
    system2(python_path, args = base_cmd, stdout = TRUE, stderr = TRUE),
    error = function(e) { paste("SYSTEM2_ERROR:", conditionMessage(e)) }
  )
  
  # Guardar log por enfermedad para depuración futura
  write_lines(as.character(py_out), glue("{yatra_dir}/yatra_run.log"))
  
  if (!file.exists(yatra_out)) {
    warning(glue("Yatra falló o no generó salida para {disease}. Se omite combinación.\nRevisa {yatra_dir}/yatra_run.log"))
    next
  }
  cat("Yatra completado para ", disease, "\n", sep = "")
  
  # === Lectura y normalización del ranking de Yatra ===
  yatra_rank <- read_tsv(yatra_out, show_col_types = FALSE)
  names_lower <- tolower(names(yatra_rank))
  col_gene  <- which(names_lower %in% c("gene","gene_name","symbol","hgnc","node"))[1]
  col_score <- which(names_lower %in% c("score","ranking","rank","prob","p"))[1]
  if (length(col_gene) == 0 || length(col_score) == 0) {
    stop("No encuentro columnas de gen/score en yatra_ranking.tsv. Columnas: ",
         paste(names(yatra_rank), collapse = ", "))
  }
  yatra_rank <- yatra_rank[, c(col_gene, col_score)]
  colnames(yatra_rank) <- c("gene","score_yatra")
  
  yatra_rank <- yatra_rank %>%
    mutate(gene = as.character(gene),
           score_yatra = suppressWarnings(as.numeric(score_yatra))) %>%
    filter(!is.na(gene)) %>%
    arrange(desc(score_yatra))
  
  # === Localizar métricas topológicas (varias rutas/convenciones) ===
  metrics_file <- first_existing(c(
    glue("output/{disease}/network/{disease}_node_metrics.tsv"),
    glue("output/{disease}/network/{disease}_node_metrics_physical.tsv"),
    glue("output/network/{disease}_node_metrics.tsv"),
    glue("output/network/{disease}_node_metrics_physical.tsv")
  ))
  
  if (length(metrics_file) > 0) {
    message("Combinando priorización Yatra + métricas topológicas")
    metrics <- read_tsv(metrics_file, show_col_types = FALSE)
    
    # Normalizar columna de gen
    nml <- tolower(names(metrics))
    if (!"gene" %in% nml) {
      idx <- which(nml %in% c("gene","gene_name","symbol","hgnc","node"))[1]
      if (length(idx) == 1) names(metrics)[idx] <- "gene"
    }
    if (!"gene" %in% tolower(names(metrics))) {
      warning("No se encontró columna de gen en métricas; se guardará solo Yatra.")
      combined <- yatra_rank %>% arrange(desc(score_yatra))
    } else {
      # Añadir columnas si faltan
      if (!"degree" %in% names(metrics))      metrics$degree <- 0
      if (!"betweenness" %in% names(metrics)) metrics$betweenness <- 0
      
      metrics <- metrics %>%
        mutate(
          degree = suppressWarnings(as.numeric(degree)),
          betweenness = suppressWarnings(as.numeric(betweenness)),
          score_topo = as.numeric(scale(replace_na(degree, 0))) +
            as.numeric(scale(replace_na(betweenness, 0)))
        ) %>%
        select(gene, score_topo)
      
      combined <- full_join(metrics, yatra_rank, by = "gene") %>%
        mutate(
          score_topo  = replace_na(score_topo, 0),
          score_yatra = replace_na(score_yatra, 0),
          score_combined = score_topo + as.numeric(scale(score_yatra))
        ) %>%
        arrange(desc(score_combined))
    }
    
    combined_out <- glue("{yatra_dir}/combined_ranking.tsv")
    write_tsv(combined, combined_out)
    cat("Ranking combinado guardado en: ", combined_out, "\n\n", sep = "")
    
  } else {
    message("No se encontraron métricas; se guarda solo el ranking de Yatra.")
    yatra_only_out <- glue("{yatra_dir}/yatra_ranking_only.tsv")
    write_tsv(yatra_rank %>% arrange(desc(score_yatra)), yatra_only_out)
  }
}
