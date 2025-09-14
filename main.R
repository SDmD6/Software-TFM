### =================== Usuario y Working Directory =================== ###
# Trabajo de TFM
# Nombre: Sofía Domenech Dauder

# Detectar automáticamente el directorio donde está este script
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  this_file <- rstudioapi::getActiveDocumentContext()$path
  if (nzchar(this_file)) {
    setwd(dirname(this_file))
  } else {
    message("Advertencia: No se pudo detectar el path del script (¿estás en consola?).")
  }
} else {
  warning("rstudioapi no está disponible. Asegúrate de estar en el directorio raíz del proyecto.")
}

### =================== PAQUETES =================== ###

# Instalar BiocManager si falta
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Paquetes Bioconductor
bioc_pkgs <- c("STRINGdb", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}

# Paquetes CRAN
cran_pkgs <- c("readr", "dplyr", "ggplot2", "ggraph", "tidygraph",
               "pheatmap", "ggvenn", "circlize", "stringr", "purrr", "digest")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

# Instala glue si falta
if (!requireNamespace("glue", quietly = TRUE)) install.packages("glue")

# Cargar librerías
suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(STRINGdb)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(pheatmap)
  library(ggvenn)
  library(circlize)
  library(digest)
  library(glue)
  library(reticulate)
})

### ======== ENTRADA DE DATOS Y FILTRADO DE ENFERMEDADES =========== ###
# Asegurar carpeta 'data/'
if (!dir.exists("data")) dir.create("data", recursive = TRUE)

# Escribir disease_pairs solo si no existe (evita sobreescritura accidental)
pairs_path <- "data/disease_pairs.tsv"
if (!file.exists(pairs_path)) {
  tibble::tibble(
    Disease = c("CMD", "CM"),
    HPO_ID  = c("HP:0006817", "HP:0003808")
  ) %>% readr::write_tsv(pairs_path)
} else {
  message("Usando data/disease_pairs.tsv existente.")
}

source("scripts/read_input.R")

### =================== ENRIQUECIMIENTO FUNCIONAL (GO + KEGG) =================== ###
source("scripts/run_enrichment.R")

### =================== REDES CON STRINGdb (funcional y física) =================== ###
source("scripts/build_network.R")

### =================== ANÁLISIS CON IGRAPH (nodos y globales) =================== ###
source("scripts/analyze_network.R")

### =================== COMPARACIONES ENTRE ENFERMEDADES Y TIPOS DE RED =================== ###
source("scripts/compare_diseases.R")

### =================== VISUALIZACIONES (GO, KEGG, redes) =================== ###
source("scripts/visual_summary.R")

### =================== PRIORIZACIÓN POR RANDOM WALK =================== ###
source("scripts/run_random_walk.R")

### =================== PRIORIZACIÓN CON YATRA =================== ###
source("scripts/prioritize_yatra.R")

### ================== CHEKEAR RESULTAODS YATRA ================== ###
source("scripts/check_yatra_results.R")

### ================== GENERAR REPORT ================== ###
source("scripts/generate_report.R")


