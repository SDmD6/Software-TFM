
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(tidyverse); library(glue)
})

write_enrichment_tsv <- function(df, path, what = "resultado") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(df) || (is.data.frame(df) && nrow(df) == 0)) {
    readr::write_tsv(tibble(mensaje = glue("Sin {what} (tabla vacía)")), path)
  } else {
    readr::write_tsv(df, path)
  }
}

map_symbols_to_entrez <- function(genes_symbol) {
  m <- bitr(genes_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  unmapped <- setdiff(unique(genes_symbol), unique(m$SYMBOL))
  list(
    entrez   = unique(m$ENTREZID),
    map_tbl  = m,
    unmapped = unmapped,
    map_rate = length(unique(m$SYMBOL)) / length(unique(genes_symbol))
  )
}

run_go_enrichment <- function(gene_list, disease) {
  mm <- map_symbols_to_entrez(gene_list)
  message(glue("[{disease}] Mapeo SYMBOL→ENTREZ: {round(100*mm$map_rate,1)}%"))
  ego <- tryCatch(
    enrichGO(
      gene          = mm$entrez,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      keyType       = "ENTREZID",
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.05,
      minGSSize     = 5,
      maxGSSize     = 500,
      readable      = TRUE
    ),
    error = function(e) { warning(glue("[{disease}] enrichGO error: {e$message}")); NULL }
  )
  out_dir_dis <- file.path("output", disease, "enrichment")
  out_disease <- file.path(out_dir_dis, paste0(disease, "_GO_enrichment.tsv"))
  out_global  <- file.path("output", "enrichment", paste0(disease, "_GO_enrichment.tsv"))
  dir.create(out_dir_dis, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(out_global), recursive = TRUE, showWarnings = FALSE)
  write_enrichment_tsv(as.data.frame(ego), out_disease, "enriquecimiento GO-BP")
  file.copy(out_disease, out_global, overwrite = TRUE)
  write_enrichment_tsv(mm$map_tbl,  file.path(out_dir_dis, paste0(disease, "_GO_mapping.tsv")),   "tabla de mapeo")
  write_enrichment_tsv(tibble(unmapped = mm$unmapped),
                       file.path(out_dir_dis, paste0(disease, "_GO_unmapped.tsv")), "genes no mapeados")
  return(ego)
}

run_kegg_enrichment <- function(gene_list, disease) {
  mm <- map_symbols_to_entrez(gene_list)
  message(glue("[{disease}] Mapeo SYMBOL→ENTREZ: {round(100*mm$map_rate,1)}% (KEGG)"))
  ekegg <- tryCatch(
    enrichKEGG(
      gene          = mm$entrez,
      organism      = "hsa",
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.05,
      minGSSize     = 5,
      maxGSSize     = 500
    ),
    error = function(e) { warning(glue("[{disease}] enrichKEGG error: {e$message}")); NULL }
  )
  if (!is.null(ekegg)) {
    ekegg <- tryCatch(setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                      error = function(e) ekegg)
  }
  out_dir_dis <- file.path("output", disease, "enrichment")
  out_disease <- file.path(out_dir_dis, paste0(disease, "_KEGG_enrichment.tsv"))
  out_global  <- file.path("output", "enrichment", paste0(disease, "_KEGG_enrichment.tsv"))
  dir.create(out_dir_dis, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(out_global), recursive = TRUE, showWarnings = FALSE)
  write_enrichment_tsv(as.data.frame(ekegg), out_disease, "enriquecimiento KEGG")
  file.copy(out_disease, out_global, overwrite = TRUE)
  write_enrichment_tsv(mm$map_tbl,  file.path(out_dir_dis, paste0(disease, "_KEGG_mapping.tsv")),   "tabla de mapeo")
  write_enrichment_tsv(tibble(unmapped = mm$unmapped),
                       file.path(out_dir_dis, paste0(disease, "_KEGG_unmapped.tsv")), "genes no mapeados")
  return(ekegg)
}

if (!exists("disease_names") || length(disease_names) < 2) {
  stop("No se encontró `disease_names` (longitud >= 2) en el entorno.")
}
if (!exists("genes_1") || !exists("genes_2")) {
  stop("Se requieren `genes_1` y `genes_2` en el entorno antes de ejecutar run_enrichment.R")
}

ego_list   <- list(); ekegg_list <- list()
gene_lists <- list(genes_1, genes_2)

for (i in seq_along(disease_names)) {
  dis <- disease_names[i]
  ego_list[[dis]]   <- run_go_enrichment(gene_lists[[i]], dis)
  ekegg_list[[dis]] <- run_kegg_enrichment(gene_lists[[i]], dis)
}

ego_1   <- ego_list[[disease_names[1]]]
ego_2   <- ego_list[[disease_names[2]]]
ekegg_1 <- ekegg_list[[disease_names[1]]]
ekegg_2 <- ekegg_list[[disease_names[2]]]

