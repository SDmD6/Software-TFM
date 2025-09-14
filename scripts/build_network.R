### =================== REDES CON STRINGdb =================== ###
library(STRINGdb)
library(tidyverse)
library(readr)

# Crear carpeta si no existe
if (!dir.exists("output/network")) dir.create("output/network", recursive = TRUE)

# Inicializar objeto STRINGdb (score mínimo 400)
string_db <- STRINGdb$new(version = "11.5", species = 9606,
                          score_threshold = 400, input_directory = "")

# Leer genes entrenados
genes_1 <- read_tsv("data/training_genes/CMD_training_genes.tsv", show_col_types = FALSE)$gene
genes_2 <- read_tsv("data/training_genes/CM_training_genes.tsv", show_col_types = FALSE)$gene

# Lista de enfermedades y genes
gene_lists <- list(
  CMD = genes_1,
  CM  = genes_2
)

# Función para construir red y guardar con distintos umbrales
build_string_network <- function(gene_list, disease_name, score_threshold, label) {
  # Mapear a IDs de STRING
  mapped <- string_db$map(data.frame(gene = gene_list), "gene", removeUnmappedRows = TRUE)
  
  # Elegir una única asignación por STRING_id (evitar duplicados en joins)
  mapped_unique <- mapped %>%
    group_by(STRING_id) %>%
    slice(1) %>%
    ungroup()
  
  # Obtener red de interacciones
  network_df <- string_db$get_interactions(mapped$STRING_id) %>%
    filter(combined_score >= score_threshold)
  
  # Convertir a nombres de genes
  network_df <- network_df %>%
    left_join(mapped_unique, by = c("from" = "STRING_id")) %>%
    rename(from_gene = gene) %>%
    left_join(mapped_unique, by = c("to" = "STRING_id")) %>%
    rename(to_gene = gene) %>%
    select(from_gene, to_gene, combined_score) %>%
    filter(!is.na(from_gene), !is.na(to_gene)) %>%
    distinct()
  
  # Guardar red
  output_file <- paste0("output/network/", disease_name, "_network_", label, ".tsv")
  write_tsv(network_df, output_file)
  cat("Red", label, "de STRING guardada en:", output_file, "\n")
}

# Ejecutar para cada enfermedad y tipo de red
for (disease in names(gene_lists)) {
  gene_list <- gene_lists[[disease]]
  
  # Red funcional (score ≥ 400)
  build_string_network(gene_list, disease, score_threshold = 400, label = "functional")
  
  # Red física (score ≥ 900)
  build_string_network(gene_list, disease, score_threshold = 900, label = "physical")
}