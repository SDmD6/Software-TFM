### =================== COMPARACIÓN ENTRE ENFERMEDADES =================== ###

library(tidyverse)
library(readr)

# Leer tabla de enfermedades
disease_table <- read_tsv("data/disease_pairs.tsv", col_names = TRUE)

# Leer métricas nodales para cada enfermedad
node_metrics <- list()
for (disease in disease_table$Disease) {
  file_path <- paste0("output/network/", disease, "_node_metrics.tsv")
  node_metrics[[disease]] <- read_tsv(file_path, show_col_types = FALSE)
}

# Priorizar genes por centralidad
prioritize_genes <- function(df, disease_name) {
  df <- df %>%
    mutate(score = as.numeric(scale(degree)) + as.numeric(scale(betweenness))) %>%
    arrange(desc(score))
  
  output_file <- paste0("output/comparison/prioritized_genes_", disease_name, ".tsv")
  write_tsv(df, output_file)
  return(df)
}

# Aplicar priorización
for (disease in disease_table$Disease) {
  node_metrics[[disease]] <- prioritize_genes(node_metrics[[disease]], disease)
}

# Leer enriquecimientos GO
enrichment_results <- list()
for (disease in disease_table$Disease) {
  file_path <- paste0("output/enrichment/", disease, "_GO_enrichment.tsv")
  enrichment_results[[disease]] <- read_tsv(file_path, show_col_types = FALSE)
}

# Comparar genes compartidos
common_genes <- intersect(node_metrics[[disease_table$Disease[1]]]$gene,
                          node_metrics[[disease_table$Disease[2]]]$gene)
cat("Genes comunes:", length(common_genes), "\n")

if (!dir.exists("output/comparison")) dir.create("output/comparison", recursive = TRUE)
write_tsv(data.frame(gene = common_genes), "output/comparison/common_genes.tsv")

# Genes exclusivos de cada enfermedad
exclusive_cmd <- setdiff(node_metrics[[disease_table$Disease[1]]]$gene,
                         node_metrics[[disease_table$Disease[2]]]$gene)
exclusive_cm  <- setdiff(node_metrics[[disease_table$Disease[2]]]$gene,
                         node_metrics[[disease_table$Disease[1]]]$gene)

write_tsv(tibble(gene = exclusive_cmd),
          file.path("output/comparison", paste0("exclusive_", disease_table$Disease[1], ".tsv")))
write_tsv(tibble(gene = exclusive_cm),
          file.path("output/comparison", paste0("exclusive_", disease_table$Disease[2], ".tsv")))


# Comparar funciones GO-BP compartidas
shared_go <- intersect(enrichment_results[[disease_table$Disease[1]]]$Description,
                       enrichment_results[[disease_table$Disease[2]]]$Description)
cat("Funciones GO-BP comunes:", length(shared_go), "\n")
write_tsv(data.frame(go_term = shared_go), "output/comparison/shared_go_terms.tsv")

# Estadísticas topológicas globales
network_summary <- function(df, disease) {
  data.frame(
    disease = disease,
    n_genes = nrow(df),
    degree_mean = mean(df$degree),
    betweenness_mean = mean(df$betweenness),
    clustering_mean = mean(df$clustering),
    n_modules = length(unique(df$module))
  )
}

summary_df <- bind_rows(
  network_summary(node_metrics[[disease_table$Disease[1]]], disease_table$Disease[1]),
  network_summary(node_metrics[[disease_table$Disease[2]]], disease_table$Disease[2])
)

write_tsv(summary_df, "output/comparison/network_summary.tsv")
print(summary_df)

# Comparar rutas KEGG compartidas
kegg_1 <- read_tsv(paste0("output/enrichment/", disease_table$Disease[1], "_KEGG_enrichment.tsv"), show_col_types = FALSE)
kegg_2 <- read_tsv(paste0("output/enrichment/", disease_table$Disease[2], "_KEGG_enrichment.tsv"), show_col_types = FALSE)

# Extraer descripciones de rutas
kegg_shared <- intersect(kegg_1$Description, kegg_2$Description)
cat("Rutas KEGG compartidas:", length(kegg_shared), "\n")

# Guardar lista de rutas compartidas
write_tsv(tibble(kegg_term = kegg_shared),
          "output/comparison/shared_kegg_paths.tsv")

# (Opcional) Limpia el alias antiguo si existe para evitar duplicados
old_alias <- "output/comparison/shared_kegg_terms.tsv"
if (file.exists(old_alias)) file.remove(old_alias)
# Estadísticas de solapamiento KEGG
kegg_overlap_summary <- tibble(
  disease_1 = disease_table$Disease[1],
  disease_2 = disease_table$Disease[2],
  n_kegg_1 = nrow(kegg_1),
  n_kegg_2 = nrow(kegg_2),
  n_shared = length(kegg_shared),
  jaccard_index = length(kegg_shared) / length(union(kegg_1$Description, kegg_2$Description))
)

write_tsv(kegg_overlap_summary, "output/comparison/kegg_overlap_summary.tsv")
print(kegg_overlap_summary)

# ---- Jaccard de rankings topológicos y Yatra (top-k) ----
suppressPackageStartupMessages({ library(readr); library(dplyr); library(glue) })

k <- 50  # ajusta aquí el k que quieras comparar

tp_path_1 <- glue("output/comparison/prioritized_genes_{disease_table$Disease[1]}.tsv")
tp_path_2 <- glue("output/comparison/prioritized_genes_{disease_table$Disease[2]}.tsv")

rows_out <- list()

# Topológico
if (file.exists(tp_path_1) && file.exists(tp_path_2)) {
  tp1 <- read_tsv(tp_path_1, show_col_types = FALSE)
  tp2 <- read_tsv(tp_path_2, show_col_types = FALSE)
  top_tp_1 <- head(tp1$gene, k); top_tp_2 <- head(tp2$gene, k)
  j_tp <- length(intersect(top_tp_1, top_tp_2)) / length(union(top_tp_1, top_tp_2))
  rows_out[[length(rows_out)+1]] <- tibble(
    method  = "topological",
    k       = k,
    disease1 = disease_table$Disease[1],
    disease2 = disease_table$Disease[2],
    jaccard = j_tp
  )
}

# Yatra (acepta dos nombres posibles de archivo)
yatra_path_1 <- glue("output/{disease_table$Disease[1]}/prioritization/yatra_ranking.tsv")
yatra_path_2 <- glue("output/{disease_table$Disease[2]}/prioritization/yatra_ranking.tsv")
if (!file.exists(yatra_path_1)) yatra_path_1 <- glue("output/{disease_table$Disease[1]}/prioritization/combined_ranking.tsv")
if (!file.exists(yatra_path_2)) yatra_path_2 <- glue("output/{disease_table$Disease[2]}/prioritization/combined_ranking.tsv")

if (file.exists(yatra_path_1) && file.exists(yatra_path_2)) {
  yr1 <- read_tsv(yatra_path_1, show_col_types = FALSE)
  yr2 <- read_tsv(yatra_path_2, show_col_types = FALSE)
  col_g1 <- if ("gene" %in% names(yr1)) "gene" else names(yr1)[1]
  col_g2 <- if ("gene" %in% names(yr2)) "gene" else names(yr2)[1]
  top_yr_1 <- head(yr1[[col_g1]], k); top_yr_2 <- head(yr2[[col_g2]], k)
  j_yr <- length(intersect(top_yr_1, top_yr_2)) / length(union(top_yr_1, top_yr_2))
  rows_out[[length(rows_out)+1]] <- tibble(
    method  = "yatra",
    k       = k,
    disease1 = disease_table$Disease[1],
    disease2 = disease_table$Disease[2],
    jaccard = j_yr
  )
}

if (length(rows_out) > 0) {
  out_tbl <- bind_rows(rows_out)
  write_tsv(out_tbl, "output/comparison/jaccard_topk_rankings.tsv")
  print(out_tbl)
} else {
  message("No se pudieron calcular Jaccard top-k (faltan rankings).")
}

### =================== RESUMEN GLOBAL CMD vs CM =================== ###
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(glue)
})

# --- Archivos de entrada ---
cmp_dir <- "output/comparison"

genes_common <- read_tsv(file.path(cmp_dir, "common_genes.tsv"), show_col_types = FALSE)
genes_cmd    <- read_tsv(file.path(cmp_dir, "exclusive_CMD.tsv"), show_col_types = FALSE)
genes_cm     <- read_tsv(file.path(cmp_dir, "exclusive_CM.tsv"), show_col_types = FALSE)

go_shared <- tryCatch(read_tsv(file.path(cmp_dir, "shared_go_terms.tsv"), show_col_types = FALSE), error = function(e) NULL)
kegg_shared <- tryCatch(read_tsv(file.path(cmp_dir, "shared_kegg_terms.tsv"), show_col_types = FALSE), error = function(e) NULL)

kegg_overlap <- read_tsv(file.path(cmp_dir, "kegg_overlap_summary.tsv"), show_col_types = FALSE)
jaccard_rank <- read_tsv(file.path(cmp_dir, "jaccard_topk_rankings.tsv"), show_col_types = FALSE)

# Leer enriquecimientos GO
go_cmd <- read_tsv("output/enrichment/CMD_GO_enrichment.tsv", show_col_types = FALSE)
go_cm  <- read_tsv("output/enrichment/CM_GO_enrichment.tsv", show_col_types = FALSE)

# Leer enriquecimientos KEGG
kegg_cmd <- read_tsv("output/enrichment/CMD_KEGG_enrichment.tsv", show_col_types = FALSE) %>% 
  filter(p.adjust <= 0.05)
kegg_cm  <- read_tsv("output/enrichment/CM_KEGG_enrichment.tsv", show_col_types = FALSE) %>% 
  filter(p.adjust <= 0.05)

# Exclusivos
exclusivos_kegg_cmd <- length(setdiff(kegg_cmd$Description, kegg_cm$Description))
exclusivos_kegg_cm  <- length(setdiff(kegg_cm$Description, kegg_cmd$Description))


# Filtrar por significancia (p.adjust ≤ 0.05)
go_cmd_sig <- go_cmd %>% filter(p.adjust <= 0.05)
go_cm_sig  <- go_cm %>% filter(p.adjust <= 0.05)

# Actualizar resumen
summary_tbl <- tibble(
  Categoria      = c("Genes", "GO-BP", "KEGG", "Rankings"),
  CMD_total      = c(nrow(genes_common) + nrow(genes_cmd),
                     nrow(go_cmd_sig),
                     nrow(kegg_cmd),
                     NA),
  CM_total       = c(nrow(genes_common) + nrow(genes_cm),
                     nrow(go_cm_sig),
                     nrow(kegg_cm),
                     NA),
  Compartidos    = c(nrow(genes_common),
                     nrow(go_shared),
                     length(intersect(kegg_cmd$Description, kegg_cm$Description)),
                     NA),
  Exclusivos_CMD = c(nrow(genes_cmd),
                     nrow(go_cmd_sig) - nrow(go_shared),
                     exclusivos_kegg_cmd,
                     NA),
  Exclusivos_CM  = c(nrow(genes_cm),
                     nrow(go_cm_sig) - nrow(go_shared),
                     exclusivos_kegg_cm,
                     NA),
  Jaccard        = c(NA,
                     round(length(unique(go_shared$go_term)) /
                             length(union(go_cmd_sig$Description, go_cm_sig$Description)), 3),
                     round(length(intersect(kegg_cmd$Description, kegg_cm$Description)) /
                             length(union(kegg_cmd$Description, kegg_cm$Description)), 3),
                     round(jaccard_rank$jaccard[1], 3))
)


# --- Guardar en TSV ---
out_file <- file.path(cmp_dir, "global_overlap_summary.tsv")
write_tsv(summary_tbl, out_file)
cat("Tabla resumen global guardada en:", out_file, "\n")

# --- Vista previa en consola ---
print(summary_tbl)
