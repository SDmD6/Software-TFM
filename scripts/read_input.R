### =================== LECTURA DE ENFERMEDADES Y GENES DESDE HPO =================== ###
library(tidyverse)
library(digest)
library(readr)

# Función auxiliar para guardar genes solo si el archivo fuente cambió
save_training_genes <- function(gene_list, disease_name, input_file) {
  output_dir <- "data/training_genes"
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  output_file <- file.path(output_dir, paste0(disease_name, "_training_genes.tsv"))
  hash_file   <- file.path(output_dir, paste0(disease_name, "_input_hash.txt"))
  
  current_hash <- digest(file = input_file, algo = "md5")
  previous_hash <- if (file.exists(hash_file)) readLines(hash_file) else ""
  
  if (!file.exists(output_file) || current_hash != previous_hash) {
    write_tsv(data.frame(gene = unique(gene_list)), output_file)
    writeLines(current_hash, hash_file)
    message("Archivo actualizado: ", output_file)
  } else {
    message("Archivo ya actualizado: ", output_file)
  }
}

# Descargar archivo HPO si no existe
hpo_file <- "data/phenotype_to_genes.txt"
if (!file.exists(hpo_file)) {
  url <- "http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt"
  download.file(url, hpo_file)
}

# Leer anotaciones HPO
hpo <- read_tsv(hpo_file, comment = "#", col_names = TRUE, show_col_types = FALSE)

# Leer lista de enfermedades y HPO IDs
disease_table <- read_tsv("data/disease_pairs.tsv", col_names = TRUE)
disease_names <- disease_table$Disease

# Inicializar lista
disease_genes <- list()

# Procesar cada enfermedad
for (i in seq_len(nrow(disease_table))) {
  disease_name <- disease_table$Disease[i]
  hpo_id <- disease_table$HPO_ID[i]
  
  genes_df <- hpo %>%
    filter(hpo_id == !!hpo_id) %>%
    select(gene_symbol) %>%
    filter(!is.na(gene_symbol)) %>%
    distinct()
  
  gene_list <- genes_df$gene_symbol
  save_training_genes(gene_list, disease_name, hpo_file)
  
  disease_genes[[disease_name]] <- gene_list
  cat("Genes ", disease_name, ": ", length(gene_list), "\n", sep = "")
}

# Variables globales necesarias para otros scripts
genes_1 <- disease_genes[[disease_names[1]]]
genes_2 <- disease_genes[[disease_names[2]]]

