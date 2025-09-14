# =================== ANALISIS TOPOLOGICO (por enfermedad) =================== #
suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(tidygraph)
  library(glue)
  library(fs)
})

# ---- helpers ----
compute_node_metrics <- function(g) {
  as_tbl_graph(g) %>%
    mutate(
      degree      = centrality_degree(mode = "all"),
      betweenness = centrality_betweenness(directed = FALSE, normalized = TRUE),
      clustering  = local_transitivity(),         # equivalente local undirected
      module      = group_walktrap()
    ) %>%
    rename(gene = name) %>%
    as_tibble()
}

compute_network_summary <- function(g, disease, network_kind) {
  comps <- igraph::components(g)
  tibble(
    disease         = disease,
    network_kind    = network_kind,                 # "physical" | "functional"
    n_nodes         = igraph::gorder(g),
    n_edges         = igraph::gsize(g),
    density         = igraph::edge_density(g, loops = FALSE),
    avg_degree      = mean(igraph::degree(g)),
    median_degree   = median(igraph::degree(g)),
    transitivity    = igraph::transitivity(g, type = "globalundirected"),
    n_components    = comps$no,
    giant_size      = ifelse(igraph::gorder(g) > 0, max(comps$csize), 0),
    giant_fraction  = ifelse(igraph::gorder(g) > 0, max(comps$csize) / igraph::gorder(g), NA_real_)
  )
}

# ---- localizar archivos de red por enfermedad ----
net_files <- dir("output", pattern = "_network_(physical|functional)\\.tsv$",
                 recursive = TRUE, full.names = TRUE)

if (length(net_files) == 0) {
  stop("No se encontraron archivos *_network_physical.tsv o *_network_functional.tsv bajo 'output/'. Ejecuta primero build_network.R.")
}

# Estructura esperada: output/<DIS>/network/<DIS>_network_<kind>.tsv
net_tbl <- tibble(path = net_files) |>
  mutate(
    file    = basename(path),
    dir     = dirname(path),
    disease = str_match(file, "^([^_]+)_network_")[,2],
    kind    = str_match(file, "_network_([^\\.]+)\\.tsv$")[,2]
  ) |>
  filter(!is.na(disease), !is.na(kind))

# ---- procesar cada red ----
by_disease <- split(net_tbl, net_tbl$disease)

for (dis in names(by_disease)) {
  rows <- by_disease[[dis]]
  outdir <- glue("output/{dis}/network")
  dir_create(outdir)
  
  summaries <- list()
  
  for (i in seq_len(nrow(rows))) {
    net_file <- rows$path[i]
    kind     <- rows$kind[i]           # physical | functional
    
    # Leer aristas
    edges <- read_tsv(net_file, show_col_types = FALSE)
    
    # Validación mínima de columnas
    needed <- c("from_gene","to_gene")
    if (!all(needed %in% colnames(edges))) {
      stop(glue("El archivo {net_file} debe contener columnas: {toString(needed)}"))
    }
    
    # Construir grafo y simplificar
    g <- igraph::graph_from_data_frame(edges[, needed], directed = FALSE)
    g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    
    # Métricas nodales
    metrics <- compute_node_metrics(g)
    
    # Guardar métricas por tipo de red
    metrics_file <- glue("{outdir}/{dis}_node_metrics_{kind}.tsv")
    write_tsv(metrics, metrics_file)
    
    # Resumen global de red
    summaries[[kind]] <- compute_network_summary(g, dis, kind)
    
    cat(glue("[OK] {dis} | {kind}: nodos={igraph::gorder(g)}; ejes={igraph::gsize(g)} -> {metrics_file}\n"))
  }
  
  # Guardar resumen por tipo y combinado
  summary_df <- bind_rows(summaries)
  write_tsv(summary_df, glue("{outdir}/{dis}_network_summary_by_kind.tsv"))
  write_tsv(summary_df, glue("{outdir}/{dis}_network_summary.tsv"))
  
  # ---- compatibilidad hacia atrás ----
  preferred <- if ("physical" %in% summary_df$network_kind) "physical" else summary_df$network_kind[1]
  src <- glue("{outdir}/{dis}_node_metrics_{preferred}.tsv")
  dst <- glue("{outdir}/{dis}_node_metrics.tsv")
  file_copy(src, dst, overwrite = TRUE)
  cat(glue("[BC] Alias compatibilidad -> {basename(dst)} -> {basename(src)}\n"))
}
