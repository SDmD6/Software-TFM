### =================== VISUALIZACIONES GENERALIZADAS =================== ###

library(tidyverse)
library(ggraph)
library(tidygraph)
library(igraph)
library(pheatmap)
library(ggvenn)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)

# --- Localizadores robustos de archivos (prioriza carpeta plana) ---
find_nodes_file <- function(disease) {
  cands <- c(
    file.path("output", "network", paste0(disease, "_node_metrics.tsv")),
    file.path("output", disease, "network", paste0(disease, "_node_metrics.tsv"))
  )
  cand <- cands[file.exists(cands)][1]
  if (is.na(cand)) stop("No encuentro métricas nodales para ", disease, " en:\n", paste(cands, collapse = "\n"))
  cand
}

find_edges_file <- function(disease, kind = c("physical","functional")) {
  kind <- match.arg(kind)
  cands <- c(
    file.path("output", "network", paste0(disease, "_network_", kind, ".tsv")),
    file.path("output", disease, "network", paste0(disease, "_network_", kind, ".tsv"))
  )
  cand <- cands[file.exists(cands)][1]
  return(cand)
}

# Elige "la mejor red disponible": física si existe; si no, funcional.
pick_best_edges <- function(disease) {
  f <- find_edges_file(disease, "physical")
  if (!is.na(f)) return(list(path = f, kind = "physical"))
  f2 <- find_edges_file(disease, "functional")
  if (!is.na(f2)) return(list(path = f2, kind = "functional"))
  stop("No encuentro red physical/functional para ", disease, " ni en output/network ni en output/<DIS>/network")
}


# ==== Utilidades comunes para enriq. (acepta enrichResult o data.frame) ====
to_df <- function(enrich_obj) {
  if (inherits(enrich_obj, c("enrichResult","gseaResult"))) as.data.frame(enrich_obj) else as.data.frame(enrich_obj)
}

# Devuelve tibble con columnas: geneID (un gen por fila), Description (término)
get_gene_pathway <- function(enrich_obj, terms_keep = NULL) {
  df <- to_df(enrich_obj) %>%
    dplyr::select(Description, geneID) %>%
    dplyr::filter(!is.na(Description), !is.na(geneID))
  if (!is.null(terms_keep)) df <- dplyr::filter(df, Description %in% terms_keep)
  df %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    dplyr::distinct(geneID, Description) %>%
    tibble::as_tibble()
}


# Leer el archivo de enfermedades
# Leer pares de enfermedades
disease_table <- readr::read_tsv("data/disease_pairs.tsv", show_col_types = FALSE)
disease_1 <- disease_table$Disease[1]
disease_2 <- disease_table$Disease[2]

# Métricas nodales (prioriza carpeta plana)
nodes_1 <- readr::read_tsv(find_nodes_file(disease_1), show_col_types = FALSE)
nodes_2 <- readr::read_tsv(find_nodes_file(disease_2), show_col_types = FALSE)

# Enriquecimiento GO (estos ya los guardas en output/enrichment)
go_1 <- readr::read_tsv(file.path("output","enrichment", paste0(disease_1, "_GO_enrichment.tsv")), show_col_types = FALSE)
go_2 <- readr::read_tsv(file.path("output","enrichment", paste0(disease_2, "_GO_enrichment.tsv")), show_col_types = FALSE)

# Resumen topológico global de la comparación
summary_df <- readr::read_tsv(file.path("output","comparison","network_summary.tsv"), show_col_types = FALSE)

# Crear carpeta si no existe
if (!dir.exists("output/plots")) dir.create("output/plots", recursive = TRUE)

### 1. Comparación de métricas topológicas globales
summary_long <- summary_df %>%
  pivot_longer(cols = -disease, names_to = "métrica", values_to = "valor") %>%
  mutate(métrica = recode(métrica,
                          degree_mean = "Grado promedio",
                          betweenness_mean = "Intermediación",
                          clustering_mean = "Agrupamiento"))

p1 <- ggplot(summary_long, aes(x = disease, y = valor, fill = disease)) +
  geom_col(position = "dodge") +
  facet_wrap(~métrica, scales = "free_y") +
  theme_minimal() +
  labs(title = "Métricas globales de red")

ggsave("output/plots/global_metrics_comparison.png", plot = p1, width = 9, height = 5)


### 2. Distribución del grado
p2 <- dplyr::bind_rows(
  dplyr::mutate(nodes_1, Enfermedad = disease_1),
  dplyr::mutate(nodes_2, Enfermedad = disease_2)
) %>%
  ggplot(aes(x = degree, fill = Enfermedad)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Distribución del grado de conectividad")

ggsave("output/plots/degree_distribution.png", plot = p2, width = 9, height = 5)


### 3. Venn de genes
venn_genes <- list()
venn_genes[[disease_1]] <- nodes_1$gene
venn_genes[[disease_2]] <- nodes_2$gene
ggvenn(venn_genes, fill_color = c("#E57373", "#64B5F6"), stroke_size = 0.5, set_name_size = 5)
ggsave("output/plots/venn_genes.png", width = 6, height = 5)

### 4. Venn de términos GO
venn_go <- list()
venn_go[[disease_1]] <- go_1$Description
venn_go[[disease_2]] <- go_2$Description
ggvenn(venn_go, fill_color = c("#E57373", "#64B5F6"), stroke_size = 0.5, set_name_size = 5)
ggsave("output/plots/venn_go_terms.png", width = 6, height = 5)

### 5. Heatmap de términos GO
suppressPackageStartupMessages({
  library(tidyverse); library(glue); library(stringr); library(fs)
})

# 1) Matriz binaria presencia/ausencia
all_terms <- union(go_1$Description, go_2$Description)
df <- tibble(
  Term = all_terms,
  !!disease_1 := as.integer(all_terms %in% go_1$Description),
  !!disease_2 := as.integer(all_terms %in% go_2$Description)
)

# 2) Long + categorías (compartidas/exclusivas)
df_long <- df %>%
  pivot_longer(-Term, names_to = "Disease", values_to = "Present")

term_cat <- df_long %>%
  group_by(Term) %>%
  summarize(
    p1 = as.integer(any(Disease == disease_1 & Present == 1)),
    p2 = as.integer(any(Disease == disease_2 & Present == 1)),
    .groups = "drop"
  ) %>%
  mutate(RowCat = case_when(
    p1 == 1 & p2 == 1 ~ "Compartidas",
    p1 == 1 & p2 == 0 ~ glue("Exclusivas {disease_1}"),
    p1 == 0 & p2 == 1 ~ glue("Exclusivas {disease_2}")
  ))

df_long <- df_long %>% left_join(term_cat, by = "Term")

# 3) (RESUMEN) Limitar nº de términos por grupo para informe
n_por_grupo <- 40  # ajusta a 20–40 para máxima legibilidad
df_long_resumen <- df_long %>%
  group_by(RowCat) %>%
  arrange(Term, .by_group = TRUE) %>%      # criterio simple y estable
  slice_head(n = n_por_grupo) %>%
  ungroup()

# 4) Orden y etiquetas envueltas (RESUMEN)
order_terms <- df_long_resumen %>%
  distinct(Term, RowCat) %>%
  arrange(RowCat, Term) %>% pull(Term)

df_long_resumen <- df_long_resumen %>%
  mutate(
    Term_wrap = str_wrap(Term, width = 60),
    Term_wrap = factor(Term_wrap, levels = str_wrap(order_terms, 60)),
    Disease   = factor(Disease, levels = c(disease_1, disease_2)),
    RowCat    = factor(RowCat, levels = c("Compartidas",
                                          glue("Exclusivas {disease_1}"),
                                          glue("Exclusivas {disease_2}")))
  )

# 5) Gráfico (objeto p_resumen reutilizable)
p_resumen <- ggplot(df_long_resumen, aes(x = Disease, y = Term_wrap, fill = factor(Present))) +
  geom_tile(color = "grey90", linewidth = 0.2, width = 0.95, height = 0.95) +
  scale_fill_manual(values = c("0" = "white", "1" = "#377eb8"),
                    labels = c("Ausente", "Presente"), name = NULL) +
  facet_grid(RowCat ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Presencia/ausencia de términos GO",
       subtitle = glue("{disease_1} vs {disease_2}"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(face = "bold"),
        axis.text.x  = element_text(face = "bold"),
        axis.text.y  = element_text(size = 7),
        legend.position = "right")

# 6) Guardado (altura proporcional) — versión RESUMEN
dir_create("output/plots")
height_cm <- max(10, nlevels(df_long_resumen$Term_wrap) * 0.20 + 6)
ggsave("output/plots/heatmap_go_terms_facetas_resumen.png", p_resumen,
       width = 26, height = height_cm, units = "cm", dpi = 300, bg = "white")


### 6. Red compartida
### ========== 6. Red compartida (con módulos nombrados) ==========

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(readr)
  library(glue)
  library(fs)
  library(stringr)
  # Para nombrar módulos con GO-BP (opcional: si no están, se usa "Hub: GEN")
  suppressWarnings({
    have_cp <- requireNamespace("clusterProfiler", quietly = TRUE) &&
      requireNamespace("org.Hs.eg.db",       quietly = TRUE)
  })
  if (have_cp) {
    library(clusterProfiler)
    library(org.Hs.eg.db)
  }
})

# ---------- Localizadores mínimos ----------
find_edges_file <- function(disease, kind = c("physical","functional")) {
  kind <- match.arg(kind)
  cands <- c(
    file.path("output", "network",          paste0(disease, "_network_", kind, ".tsv")),
    file.path("output", disease, "network", paste0(disease, "_network_", kind, ".tsv"))
  )
  cand <- cands[file.exists(cands)][1]
  ifelse(length(cand) == 0, NA_character_, cand)
}
pick_best_edges <- function(disease) {
  f <- find_edges_file(disease, "physical")
  if (!is.na(f)) return(list(path = f, kind = "physical"))
  f2 <- find_edges_file(disease, "functional")
  if (!is.na(f2)) return(list(path = f2, kind = "functional"))
  stop("No encuentro red physical/functional para ", disease)
}

# ---------- Helper: nombres legibles para módulos ----------
quiet_try <- function(expr) suppressWarnings(suppressMessages(try(expr, silent = TRUE)))

# Criterio: término GO-BP más significativo; si no hay, "Hub: <GEN>"
make_module_labels <- function(nodes_tbl, wrap = 22) {
  stopifnot(all(c("name","module","degree") %in% names(nodes_tbl)))
  
  # niveles en orden estable
  levs <- levels(as.factor(nodes_tbl$module))
  mods <- levs
  raw_labels <- character(length(mods))
  
  for (i in seq_along(mods)) {
    idx   <- which(nodes_tbl$module == mods[i])
    genes <- nodes_tbl$name[idx]
    hub   <- genes[ which.max(nodes_tbl$degree[idx]) ]
    label <- paste0("Hub: ", hub)
    
    if (isTRUE(have_cp) && length(genes) >= 3) {
      conv <- quiet_try(clusterProfiler::bitr(
        genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db
      ))
      if (!inherits(conv, "try-error") && !is.null(conv) && nrow(conv) > 0) {
        ego <- quiet_try(clusterProfiler::enrichGO(
          gene = unique(conv$ENTREZID), OrgDb = org.Hs.eg.db, ont = "BP",
          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
          readable = TRUE
        ))
        if (!inherits(ego, "try-error") && !is.null(ego)) {
          edf <- as.data.frame(ego)
          if (!is.null(edf) && nrow(edf) > 0 && "Description" %in% names(edf)) {
            label <- edf$Description[1]
          }
        }
      }
    }
    raw_labels[i] <- str_wrap(label, width = wrap)
  }
  
  # Desambiguar duplicadas: añadir id de módulo y asegurar unicidad
  dup <- duplicated(raw_labels) | duplicated(raw_labels, fromLast = TRUE)
  lab <- raw_labels
  if (any(dup)) {
    lab[dup] <- paste0(raw_labels[dup], " [M", mods[dup], "]")
    lab <- make.unique(lab, sep = " · ")
  }
  
  tibble(
    module       = factor(mods, levels = levs),
    module_label = factor(lab,  levels = lab)  # niveles únicos -> leyenda estable
  )
}

# ---------- Entradas: enfermedades y genes compartidos ----------
pairs <- readr::read_tsv("data/disease_pairs.tsv", show_col_types = FALSE)
disease_1 <- pairs$Disease[1]
disease_2 <- pairs$Disease[2]

if (!file.exists("output/comparison/common_genes.tsv")) {
  stop("Falta output/comparison/common_genes.tsv. Ejecuta compare_diseases.R antes.")
}
common_genes <- readr::read_tsv("output/comparison/common_genes.tsv",
                                show_col_types = FALSE)$gene

# ---------- Red base: mejor red de disease_1 (physical > functional) ----------
best1   <- pick_best_edges(disease_1)
edges_1 <- readr::read_tsv(best1$path, show_col_types = FALSE)
stopifnot(all(c("from_gene","to_gene") %in% names(edges_1)))

# Subred inducida por genes compartidos
edges_shared <- edges_1 %>%
  dplyr::filter(from_gene %in% common_genes, to_gene %in% common_genes) %>%
  dplyr::distinct(from_gene, to_gene)

if (nrow(edges_shared) == 0L) {
  stop("No hay aristas entre los genes compartidos en la red seleccionada (", best1$kind, ").")
}

# ---------- Grafo, módulos, degree ----------
g  <- igraph::graph_from_data_frame(edges_shared, directed = FALSE) |> igraph::simplify()
tg <- tidygraph::as_tbl_graph(g) |>
  dplyr::mutate(module = as.factor(group_walktrap()),
                degree  = centrality_degree())

# ---------- Nombres “bonitos” de módulo ----------
lab_df   <- make_module_labels(tidygraph::as_tibble(tg), wrap = 22)
tg_named <- tg %>% dplyr::left_join(lab_df, by = "module")

# ---------- Figura ----------
p_shared <- ggraph::ggraph(tg_named, layout = "fr") +
  ggraph::geom_edge_link(alpha = 0.25, colour = "grey70") +
  ggraph::geom_node_point(aes(size = degree, color = module_label)) +
  ggraph::geom_node_text(aes(label = name), family = "sans", repel = TRUE, size = 3) +
  scale_size_continuous(range = c(2, 6), name = "Grado") +
  guides(color = guide_legend(title = "Módulo (nombre)", override.aes = list(size = 5)),
         size  = guide_legend(order = 2)) +
  ggraph::theme_graph(base_size = 11, base_family = "sans") +
  labs(title    = glue::glue("Red de genes compartidos: {disease_1} ∩ {disease_2}"),
       subtitle = glue::glue("Usando red {best1$kind} de {disease_1}"))

# ---------- Guardado ----------
fs::dir_create("output/plots")
if (capabilities("cairo")) {
  ggplot2::ggsave("output/plots/shared_network_modules.png", p_shared,
                  device = "png", type = "cairo",
                  width = 10, height = 8, dpi = 300, bg = "white")
} else {
  ggplot2::ggsave("output/plots/shared_network_modules.png", p_shared,
                  width = 10, height = 8, dpi = 300, bg = "white")
}


# =========7.Red anotada por enfermedad (con nombres de nodos) =========

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(readr)
  library(glue)
  library(fs)
  library(clusterProfiler)   # si no está, el script hace fallback a "Hub: GEN"
  library(org.Hs.eg.db)
})

# ---------- Localizadores de archivos ----------
find_edges_file <- function(disease, kind = c("physical","functional")) {
  kind <- match.arg(kind)
  cands <- c(
    file.path("output", "network",             paste0(disease, "_network_", kind, ".tsv")),
    file.path("output", disease, "network",    paste0(disease, "_network_", kind, ".tsv"))
  )
  cand <- cands[file.exists(cands)][1]
  ifelse(length(cand) == 0, NA_character_, cand)
}

# Prioriza physical; si no existe, cae a functional
pick_best_edges <- function(disease) {
  f <- find_edges_file(disease, "physical")
  if (!is.na(f)) return(list(path = f, kind = "physical"))
  f2 <- find_edges_file(disease, "functional")
  if (!is.na(f2)) return(list(path = f2, kind = "functional"))
  stop("No encuentro red physical/functional para ", disease,
       " ni en output/network ni en output/<ENF>/network")
}

# ---------- Helper: pon nombre legible a cada módulo ----------
# Criterio: GO-BP más significativo; si no hay, "Hub: <GEN>"
make_module_labels <- function(nodes_tbl, wrap = 22) {
  stopifnot(all(c("name","module","degree") %in% names(nodes_tbl)))
  mods <- levels(nodes_tbl$module)
  out  <- character(length(mods))

  have_cp <- requireNamespace("clusterProfiler", quietly = TRUE) &&
             requireNamespace("org.Hs.eg.db", quietly = TRUE)

  for (i in seq_along(mods)) {
    idx   <- which(nodes_tbl$module == mods[i])
    genes <- nodes_tbl$name[idx]
    # Fallback por defecto
    hub   <- genes[which.max(nodes_tbl$degree[idx])]
    label <- paste0("Hub: ", hub)

    if (have_cp && length(genes) >= 3) {
      conv <- try(clusterProfiler::bitr(genes, fromType = "SYMBOL",
                                        toType   = "ENTREZID",
                                        OrgDb    = org.Hs.eg.db),
                  silent = TRUE)
      if (!inherits(conv, "try-error") && !is.null(conv) && nrow(conv) > 0) {
        ego <- try(clusterProfiler::enrichGO(
          gene          = unique(conv$ENTREZID),
          OrgDb         = org.Hs.eg.db,
          ont           = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff  = 0.05,
          qvalueCutoff  = 0.2,
          readable      = TRUE
        ), silent = TRUE)
        if (!inherits(ego, "try-error") && !is.null(ego)) {
          edf <- as.data.frame(ego)
          if (!is.null(edf) && nrow(edf) > 0 && "Description" %in% names(edf)) {
            label <- edf$Description[1]
          }
        }
      }
    }
    out[i] <- stringr::str_wrap(label, width = wrap)
  }

  tibble(
    module       = factor(mods, levels = mods),
    module_label = factor(out,  levels = out)  # fija el mapeo color↔etiqueta
  )
}

# ---------- Figura anotada por enfermedad ----------
plot_disease_network_annotated <- function(disease,
                                           kind = c("auto","physical","functional"),
                                           only_giant = TRUE,
                                           layout = "fr",
                                           out_dir = file.path("output", disease, "plots")) {
  kind <- match.arg(kind)
  fs::dir_create(out_dir)

  # 1) Localizar red
  if (kind == "auto") {
    best <- pick_best_edges(disease)
    edges_path <- best$path
    kind_used  <- best$kind
  } else {
    f <- find_edges_file(disease, kind)
    if (is.na(f)) stop("No encuentro red ", kind, " para ", disease)
    edges_path <- f
    kind_used  <- kind
  }

  # 2) Grafo (opcional: solo componente gigante para legibilidad)
  edges <- readr::read_tsv(edges_path, show_col_types = FALSE)
  stopifnot(all(c("from_gene","to_gene") %in% names(edges)))
  g <- igraph::graph_from_data_frame(edges[, c("from_gene","to_gene")], directed = FALSE) |>
       igraph::simplify()
  if (only_giant && igraph::gorder(g) > 0) {
    comps <- igraph::components(g)
    g <- igraph::induced_subgraph(g, which(comps$membership == which.max(comps$csize)))
  }
  if (igraph::gorder(g) == 0) stop("La red de ", disease, " está vacía tras el filtrado.")

  # 3) Métricas + módulos y nombres de módulo
  tg <- tidygraph::as_tbl_graph(g) |>
    dplyr::mutate(
      degree = centrality_degree(mode = "all"),
      module = as.factor(group_louvain())
    )
  lab_df <- make_module_labels(tidygraph::as_tibble(tg), wrap = 22)
  tg <- tg |> dplyr::left_join(lab_df, by = "module")

  # 4) Figura
  p <- ggraph::ggraph(tg, layout = layout) +
    ggraph::geom_edge_link(alpha = 0.25, colour = "grey60") +
    ggraph::geom_node_point(aes(size = degree, color = module_label)) +
    ggraph::geom_node_text(aes(label = name), family = "sans", repel = TRUE, size = 3) +
    scale_size_continuous(range = c(2, 6), name = "Grado") +
    guides(color = guide_legend(title = "Módulo (nombre)", override.aes = list(size = 5)),
           size  = guide_legend(order = 2)) +
    ggraph::theme_graph(base_size = 12, base_family = "sans") +
    labs(
      title    = glue("Red anotada: {disease}"),
      subtitle = glue("Usando red {kind_used} de {disease}"),
      x = NULL, y = NULL
    )

  # 5) Guardado (Cairo evita problemas de fuentes en Windows)
  pngf <- file.path(out_dir, paste0(disease, "_network_annotated.png"))
  pdff <- file.path(out_dir, paste0(disease, "_network_annotated.pdf"))
  use_cairo <- capabilities("cairo")

  if (use_cairo) {
    ggplot2::ggsave(pngf, p, device = "png", type = "cairo",
                    width = 11, height = 8, dpi = 300, bg = "white")
    ggplot2::ggsave(pdff, p, device = grDevices::cairo_pdf,
                    width = 11, height = 8)
  } else {
    ggplot2::ggsave(pngf, p, width = 11, height = 8, dpi = 300, bg = "white")
    message("Cairo no disponible: omito el PDF para evitar problemas de fuentes.")
  }
  message(glue("Figura guardada: {pngf}{if (use_cairo) paste0('  |  ', pdff) else ''}"))
}

# ---------- Ejecutar para el par de enfermedades ----------
pairs <- readr::read_tsv("data/disease_pairs.tsv", show_col_types = FALSE)
disease_1 <- pairs$Disease[1]
disease_2 <- pairs$Disease[2]

for (d in c(disease_1, disease_2)) {
  plot_disease_network_annotated(d, kind = "auto", only_giant = TRUE)
}

### =================== COMPARACIÓN FUNCIONAL VS FÍSICA POR ENFERMEDAD =================== ###
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(ggvenn); library(fs)
})

# Cargar tabla de enfermedades si no existe en el entorno
if (!exists("disease_table")) {
  disease_table <- readr::read_tsv("data/disease_pairs.tsv", show_col_types = FALSE)
}

# helper: busca archivo primero en subcarpeta por enfermedad, luego en carpeta plana output/network
find_netfile <- function(disease, stem) {
  cands <- c(
    file.path("output", disease, "network", paste0(disease, "_", stem, ".tsv")),
    file.path("output", "network",            paste0(disease, "_", stem, ".tsv"))
  )
  cand <- cands[file.exists(cands)][1]
  if (is.na(cand)) return(NA_character_)
  cand
}

for (disease in disease_table$Disease) {
  net_dir <- file.path("output", disease, "network")
  fs::dir_create(net_dir)
  
  # --- nodos funcional y físico (métricas) ---
  f_fun <- find_netfile(disease, "node_metrics_functional")
  f_phy <- find_netfile(disease, "node_metrics_physical")   # <-- aquí estaba el error
  
  if (is.na(f_fun) || is.na(f_phy)) {
    warning("Faltan métricas nodales functional/physical para ", disease,
            ". fun=", f_fun, " phy=", f_phy, ". Se omite esta enfermedad.")
    next
  }
  
  nodes_fun <- readr::read_tsv(f_fun, show_col_types = FALSE)
  nodes_phy <- readr::read_tsv(f_phy, show_col_types = FALSE)
  
  # === VENN DE GENES (funcional vs física) ===
  if (!("gene" %in% names(nodes_fun) && "gene" %in% names(nodes_phy))) {
    warning("No encuentro columna 'gene' en métricas de ", disease, ". Venn omitido.")
  } else {
    venn_genes_fp <- list(Functional = nodes_fun$gene, Physical = nodes_phy$gene)
    p_venn_fp <- ggvenn(venn_genes_fp,
                        fill_color = c("#FFA07A", "#20B2AA"),
                        stroke_size = 0.5, set_name_size = 4)
    ggsave(filename = file.path(net_dir, paste0(disease, "_venn_functional_vs_physical.png")),
           plot = p_venn_fp, width = 5, height = 5, dpi = 300, bg = "white")
  }
  
  # === BARPLOT DE MÉTRICAS GLOBALES (robusto a nombres) ===
  f_sum <- find_netfile(disease, "network_summary_by_kind")
  if (is.na(f_sum)) {
    warning("No existe resumen por tipo para ", disease, ". Barplot omitido.")
    next
  }
  
  summary_df <- readr::read_tsv(f_sum, show_col_types = FALSE) %>%
    dplyr::rename(
      type                   = dplyr::any_of("network_kind"),
      largest_component_size = dplyr::any_of("giant_size")
    ) %>%
    { if (!"diameter" %in% names(.)) dplyr::mutate(., diameter = NA_real_) else . } %>%
    { if ("type" %in% names(.)) dplyr::mutate(., type = factor(type, levels = c("functional","physical"))) else . }
  
  metricas_posibles  <- c("n_nodes","n_edges","avg_degree","density","diameter",
                          "n_components","largest_component_size")
  metricas_presentes <- intersect(metricas_posibles, names(summary_df))
  
  if (length(metricas_presentes) == 0L || !"type" %in% names(summary_df)) {
    warning("Estructura inesperada en resumen de ", disease, ". Barplot omitido.")
    next
  }
  
  metrics_long <- summary_df %>%
    tidyr::pivot_longer(cols = dplyr::all_of(metricas_presentes),
                        names_to = "métrica", values_to = "valor") %>%
    dplyr::mutate(métrica = dplyr::recode(métrica,
                                          n_nodes = "Nodos",
                                          n_edges = "Aristas",
                                          avg_degree = "Grado medio",
                                          density = "Densidad",
                                          diameter = "Diámetro",
                                          n_components = "Nº Componentes",
                                          largest_component_size = "Componente principal"
    ))
  
  p_metrics <- ggplot2::ggplot(metrics_long, ggplot2::aes(x = `métrica`, y = valor, fill = type)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("Red funcional vs física:", disease),
                  x = NULL, y = NULL, fill = "Tipo de red") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::scale_fill_manual(values = c("functional" = "#FFA07A", "physical" = "#20B2AA"))
  
  ggplot2::ggsave(filename = file.path(net_dir, paste0(disease, "_barplot_functional_vs_physical.png")),
                  plot = p_metrics, width = 8, height = 4, dpi = 300, bg = "white")
}  

### =================== VISUALIZACIONES DE ENRIQUECIMIENTO KEGG =================== ###

# Leer enriquecimiento KEGG (pueden ser data.frames exportados de clusterProfiler::as.data.frame)
kegg_1 <- read_tsv(paste0("output/enrichment/", disease_1, "_KEGG_enrichment.tsv"), show_col_types = FALSE)
kegg_2 <- read_tsv(paste0("output/enrichment/", disease_2, "_KEGG_enrichment.tsv"), show_col_types = FALSE)

# === Barplot KEGG por enfermedad (Top 10 por p.adjust) ===
top_kegg_1 <- kegg_1 %>% dplyr::slice_min(p.adjust, n = 10) %>% dplyr::mutate(Disease = disease_1)
top_kegg_2 <- kegg_2 %>% dplyr::slice_min(p.adjust, n = 10) %>% dplyr::mutate(Disease = disease_2)
kegg_combined <- dplyr::bind_rows(top_kegg_1, top_kegg_2)

p_kegg_bar <- ggplot(kegg_combined,
                     aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = Disease)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 rutas KEGG enriquecidas por enfermedad", x = NULL, y = "-log10(p.adjust)")

ggsave("output/plots/kegg_barplot.png", plot = p_kegg_bar, width = 9, height = 5, dpi = 300)

# === Venn KEGG ===
venn_kegg <- list()
venn_kegg[[disease_1]] <- kegg_1$Description
venn_kegg[[disease_2]] <- kegg_2$Description
p_venn_kegg <- ggvenn(venn_kegg, fill_color = c("#E57373", "#64B5F6"), stroke_size = 0.5, set_name_size = 5)
ggsave("output/plots/venn_kegg_terms.png", plot = p_venn_kegg, width = 6, height = 5, dpi = 300)

# === Heatmap KEGG (presencia/ausencia de rutas) ===
# --- Construcción matriz binaria ---
all_kegg <- union(kegg_1$Description, kegg_2$Description)
df <- tibble(
  Term = all_kegg,
  !!disease_1 := as.integer(all_kegg %in% kegg_1$Description),
  !!disease_2 := as.integer(all_kegg %in% kegg_2$Description)
)

# --- Long + categorías de fila (compartidas/exclusivas) ---
df_long <- df %>%
  pivot_longer(-Term, names_to = "Disease", values_to = "Present")

term_cat <- df_long %>%
  group_by(Term) %>%
  summarize(
    p1 = as.integer(any(Disease == disease_1 & Present == 1)),
    p2 = as.integer(any(Disease == disease_2 & Present == 1)),
    .groups = "drop"
  ) %>%
  mutate(RowCat = case_when(
    p1 == 1 & p2 == 1 ~ "Compartidas",
    p1 == 1 & p2 == 0 ~ glue("Exclusivas {disease_1}"),
    p1 == 0 & p2 == 1 ~ glue("Exclusivas {disease_2}")
  )) %>% select(Term, RowCat)

df_long <- df_long %>% left_join(term_cat, by = "Term")

# --- (Opcional) limitar nº de términos por grupo para legibilidad ---
max_terms_per_group <- 35  # AJUSTA: 25–50 suele ir bien
if (!is.null(max_terms_per_group)) {
  keep_terms <- df_long %>%
    distinct(Term, RowCat) %>%
    group_by(RowCat) %>%
    arrange(Term, .by_group = TRUE) %>%       # aquí puedes cambiar el criterio
    slice_head(n = max_terms_per_group) %>%
    ungroup() %>% pull(Term)
  df_long <- df_long %>% filter(Term %in% keep_terms)
}

# --- Orden: por grupo y alfabético dentro de grupo ---
order_terms <- df_long %>%
  distinct(Term, RowCat) %>%
  mutate(RowCat = factor(RowCat,
                         levels = c("Compartidas",
                                    glue("Exclusivas {disease_1}"),
                                    glue("Exclusivas {disease_2}")))) %>%
  arrange(RowCat, Term) %>% pull(Term)

df_long <- df_long %>%
  mutate(
    Term_wrap = str_wrap(Term, width = 55),
    Term_wrap = factor(Term_wrap, levels = str_wrap(order_terms, 55)),
    Disease = factor(Disease, levels = c(disease_1, disease_2)),
    RowCat = factor(RowCat, levels = c("Compartidas",
                                       glue("Exclusivas {disease_1}"),
                                       glue("Exclusivas {disease_2}")))
  )

# --- Gráfico ---
p <- ggplot(df_long, aes(x = Disease, y = Term_wrap, fill = factor(Present))) +
  geom_tile(color = "grey85", linewidth = 0.3, width = 0.95, height = 0.95) +
  scale_fill_manual(
    values = c("0" = "white", "1" = "#F28E2B"),
    labels = c("Ausente", "Presente"), name = NULL
  ) +
  facet_grid(RowCat ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Presencia/ausencia de rutas KEGG",
    subtitle = glue("{disease_1} vs {disease_2}: rutas compartidas y exclusivas"),
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.text.y = element_text(face = "bold"),
    axis.text.x  = element_text(face = "bold"),
    axis.text.y  = element_text(size = 8),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = 4)),
    plot.subtitle = element_text(color = "grey30"),
    plot.margin = margin(t = 6, r = 8, b = 6, l = 6)
  ) +
  coord_cartesian(clip = "off")

# --- Exportación (apaisado) ---
dir_create("output/plots")
# Altura proporcional al nº de términos mostrados
n_terms <- nlevels(df_long$Term_wrap)
height_cm <- max(10, n_terms * 0.22 + 6)  # escala suave y piso mínimo
ggsave("output/plots/heatmap_kegg_terms_facetas.png", p,
       width = 28, height = height_cm, units = "cm", dpi = 300, bg = "white")

# === Chord plot KEGG (genes ↔ rutas compartidas) ===

kegg_shared <- intersect(kegg_1$Description, kegg_2$Description)

if (length(kegg_shared) == 0) {
  message("No hay rutas KEGG compartidas entre ", disease_1, " y ", disease_2, ". No se genera chord plot.")
} else {
  # Construir tablas gene↔pathway para ambas enfermedades, restringidas a las compartidas
  kegg_map_1 <- get_gene_pathway(kegg_1, kegg_shared) %>% dplyr::mutate(Disease = disease_1)
  kegg_map_2 <- get_gene_pathway(kegg_2, kegg_shared) %>% dplyr::mutate(Disease = disease_2)
  kegg_map   <- dplyr::bind_rows(kegg_map_1, kegg_map_2) %>% dplyr::distinct()
  
  # (Opcional) limitar a términos más informativos (Top 10 por nº de genes conectados) para legibilidad
  term_counts <- kegg_map %>% dplyr::count(Description, sort = TRUE)
  n_terms <- min(10L, nrow(term_counts))
  top_terms <- term_counts %>% dplyr::slice_head(n = n_terms) %>% dplyr::pull(Description)
  kegg_map_top <- kegg_map %>% dplyr::filter(Description %in% top_terms)
  
  # Construir matriz (genes x términos) como tabla de contingencia:
  mat <- xtabs(~ geneID + Description, data = kegg_map_top)  # 0/1
  
  # Dibujar chord (cada 1 se convierte en un enlace gene ↔ término)
  # circlize acepta directamente un data.frame de pares (from,to) o una matriz:
  png("output/plots/chordplot_kegg_shared.png", width = 2000, height = 2000, res = 260)
  circos.clear()
  # Espaciado para legibilidad
  all_from <- rownames(mat); all_to <- colnames(mat)
  circos.par(gap.after = c(rep(1, length(all_from)-1), 8, rep(1, length(all_to)-1), 8))
  chordDiagram(mat, transparency = 0.25, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.06))
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.name,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.45)
  }, bg.border = NA)
  title(main = sprintf("Genes ↔ rutas KEGG compartidas (Top %d términos)", length(top_terms)))
  dev.off()
}


# ========= Figura: Subred top-N priorizados por Yatra con módulos (leyenda) =========
suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(glue)
  library(scales)
})

# Directorios de salida
dir.create("output/comparison/plots", recursive = TRUE, showWarnings = FALSE)

draw_top_yatra_subgraph <- function(
    disease,
    topN = 20,
    # Cambia a la red funcional si lo prefieres:
    net_path = glue("output/network/{disease}_network_physical.tsv"),
    out_dir  = glue("output/{disease}/plots"),
    expand_if_sparse = TRUE  # si hay pocas aristas internas, amplía a 1-vecino
) {
  if (!file.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # 1) Ranking Yatra (prioridad al combinado; si no existe, usa yatra_ranking)
  rk_candidates <- c(
    glue("output/{disease}/prioritization/combined_ranking.tsv"),
    glue("output/{disease}/prioritization/yatra_ranking.tsv")
  )
  rk_file <- rk_candidates[file.exists(rk_candidates)][1]
  if (is.na(rk_file)) {
    warning(glue("No hay ranking Yatra para {disease}."))
    return(invisible(NULL))
  }
  
  rk <- readr::read_tsv(rk_file, show_col_types = FALSE)
  
  # Normaliza nombres (acepta 'gene' o 'gene_name'; 'score_yatra' o 'score')
  nml <- tolower(names(rk))
  if (!"gene" %in% nml) {
    idxg <- which(nml %in% c("gene","gene_name","symbol","node"))[1]
    if (!is.na(idxg)) names(rk)[idxg] <- "gene"
  }
  if (!"score_yatra" %in% names(rk)) {
    if ("score" %in% names(rk)) rk <- rk %>% dplyr::rename(score_yatra = score)
  }
  stopifnot(all(c("gene","score_yatra") %in% names(rk)))
  
  rk_top <- rk %>%
    arrange(desc(score_yatra)) %>%
    slice_head(n = topN)
  
  # 2) Red (física por defecto)
  if (!file.exists(net_path)) {
    stop(glue("No existe la red para {disease}: {net_path}"))
  }
  net_df <- readr::read_tsv(net_path, show_col_types = FALSE)
  
  # Normaliza columnas de aristas a (from, to, score)
  col_from  <- intersect(c("from_gene","from","source","u","Source"), names(net_df))[1]
  col_to    <- intersect(c("to_gene","to","target","v","Target"), names(net_df))[1]
  col_score <- intersect(c("combined_score","score","weight","w","Weight"), names(net_df))[1]
  if (any(is.na(c(col_from, col_to, col_score)))) {
    stop("No encuentro columnas equivalentes a (from, to, score) en la red.")
  }
  
  edges <- net_df %>%
    transmute(from = .data[[col_from]],
              to   = .data[[col_to]],
              w    = as.numeric(.data[[col_score]])) %>%
    filter(!is.na(from), !is.na(to), from != to)
  
  genes_top <- rk_top$gene
  
  # 3) Subred inducida interna; si queda demasiado pequeña, amplía a 1-vecino
  sub_edges <- edges %>% filter(from %in% genes_top & to %in% genes_top)
  expanded <- FALSE
  if (nrow(sub_edges) < 3 && expand_if_sparse) {
    sub_edges <- edges %>% filter(from %in% genes_top | to %in% genes_top)
    expanded <- TRUE
  }
  
  g <- graph_from_data_frame(sub_edges, directed = FALSE)
  
  if (vcount(g) == 0) {
    warning(glue("Subred vacía para {disease} (top {topN})."))
    return(invisible(NULL))
  }
  
  # Alinea el orden de nodos y asigna atributos
  V(g)$name <- as.character(V(g)$name)
  
  # Atributo: score_yatra (0 si no está en topN)
  score_map <- rk_top %>% select(gene, score_yatra)
  V(g)$score_yatra <- score_map$score_yatra[match(V(g)$name, score_map$gene)]
  V(g)$score_yatra[is.na(V(g)$score_yatra)] <- 0
  
  # 4) Comunidades (módulos) con Louvain; fallback a componentes si falla
  module <- tryCatch(
    cluster_louvain(g)$membership,
    error = function(e) components(g)$membership
  )
  V(g)$module <- factor(module)
  
  # 5) Tamaño de nodo por score (escala suave), grosor de arista por peso
  node_size <- rescale(V(g)$score_yatra, to = c(3, 10), from = range(V(g)$score_yatra, na.rm = TRUE))
  node_size[is.na(node_size)] <- 3
  
  ew <- edge_attr(g, "w")
  edge_w <- rescale(ew, to = c(0.3, 2), from = range(ew, na.rm = TRUE))
  
  # 6) Título y subtítulo
  subt <- if (expanded)
    "Ampliado a 1-vecino por baja conectividad interna"
  else
    "Subred inducida interna (solo aristas entre genes top-N)"
  
  # 7) Plot con leyenda de módulos
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = edge_w), alpha = 0.35, colour = "grey40") +
    geom_node_point(aes(color = module, size = node_size), alpha = 0.95) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_size_identity() +                # el tamaño ya está en unidades
    scale_edge_width_identity() +          # el grosor ya está en unidades
    guides(color = guide_legend(title = "Módulo (Louvain)")) +
    labs(
      title = glue("{disease} · Subred top-{topN} priorizados por Yatra"),
      subtitle = paste(subt, "· Tamaño nodo ~ score Yatra · Grosor arista ~ peso STRING"),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid = element_blank()
    )
  
  # 8) Guardar PNG y PDF
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  pngf <- glue("{out_dir}/{disease}_yatra_top{topN}_subnetwork.png")
  pdff <- glue("{out_dir}/{disease}_yatra_top{topN}_subnetwork.pdf")
  ggsave(pngf, p, width = 9, height = 6, dpi = 200)
  ggsave(pdff, p, width = 9, height = 6, device = "pdf")
  
  message(glue("Subred Yatra top-{topN} para {disease} guardada en:\n- {pngf}\n- {pdff}"))
}

# Ejecutar para tus enfermedades (ajusta si lees dinámicamente desde disease_pairs.tsv)
for (d in c("CMD","CM")) {
  draw_top_yatra_subgraph(d, topN = 20)  # cambia topN si lo necesitas
}

# ========= Tabla: Resumen por módulo (top-N Yatra) con GO/KEGG =========
suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(glue)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

summarize_top_yatra_modules <- function(
    disease,
    topN = 20,
    net_path = glue("output/network/{disease}_network_physical.tsv"),
    ranking_candidates = c(
      glue("output/{disease}/prioritization/combined_ranking.tsv"),
      glue("output/{disease}/prioritization/yatra_ranking.tsv")
    ),
    expand_if_sparse = TRUE,     # igual que en la figura
    k_top_terms = 3              # nº de términos GO/KEGG a resumir por módulo
){
  # ------ 0) Cargar ranking ------
  rk_file <- ranking_candidates[file.exists(ranking_candidates)][1]
  if (is.na(rk_file)) {
    warning(glue("No hay ranking Yatra para {disease}. No se puede resumir módulos."))
    return(invisible(NULL))
  }
  rk <- readr::read_tsv(rk_file, show_col_types = FALSE)
  nml <- tolower(names(rk))
  if (!"gene" %in% nml) {
    idxg <- which(nml %in% c("gene","gene_name","symbol","node"))[1]
    if (!is.na(idxg)) names(rk)[idxg] <- "gene"
  }
  if (!"score_yatra" %in% names(rk)) {
    if ("score" %in% names(rk)) rk <- rk %>% dplyr::rename(score_yatra = score)
  }
  stopifnot(all(c("gene","score_yatra") %in% names(rk)))
  
  rk_top <- rk %>% arrange(desc(score_yatra)) %>% slice_head(n = topN)
  
  # ------ 1) Cargar red y obtener subred ------
  if (!file.exists(net_path)) stop(glue("No existe la red para {disease}: {net_path}"))
  net_df <- readr::read_tsv(net_path, show_col_types = FALSE)
  
  col_from  <- intersect(c("from_gene","from","source","u","Source"), names(net_df))[1]
  col_to    <- intersect(c("to_gene","to","target","v","Target"), names(net_df))[1]
  col_score <- intersect(c("combined_score","score","weight","w","Weight"), names(net_df))[1]
  if (any(is.na(c(col_from, col_to, col_score)))) stop("No encuentro columnas (from,to,score) en la red.")
  
  edges <- net_df %>%
    transmute(from = .data[[col_from]], to = .data[[col_to]],
              w = as.numeric(.data[[col_score]])) %>%
    filter(!is.na(from), !is.na(to), from != to)
  
  genes_top <- rk_top$gene
  sub_edges <- edges %>% filter(from %in% genes_top & to %in% genes_top)
  expanded <- FALSE
  if (nrow(sub_edges) < 3 && expand_if_sparse) {
    sub_edges <- edges %>% filter(from %in% genes_top | to %in% genes_top)
    expanded <- TRUE
  }
  
  g <- graph_from_data_frame(sub_edges, directed = FALSE)
  if (vcount(g) == 0) {
    warning(glue("Subred vacía para {disease} (top {topN})."))
    return(invisible(NULL))
  }
  V(g)$name <- as.character(V(g)$name)
  
  # ------ 2) Comunidades (Louvain) ------
  module <- tryCatch(
    cluster_louvain(g)$membership,
    error = function(e) components(g)$membership
  )
  V(g)$module <- as.integer(module)
  
  # ------ 3) Preparar listas por módulo ------
  df_nodes <- tibble(
    gene   = V(g)$name,
    module = V(g)$module
  ) %>%
    left_join(rk_top %>% select(gene, score_yatra), by = "gene") %>%
    mutate(score_yatra = replace_na(score_yatra, 0))
  
  # guardar tabla de genes por módulo (para trazabilidad)
  genes_out <- glue("output/{disease}/prioritization/module_genes_top{topN}.tsv")
  dir.create(dirname(genes_out), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(df_nodes %>% arrange(module, desc(score_yatra)), genes_out)
  
  # ------ 4) Enriquecimiento por módulo (GO-BP y KEGG) ------
  # Mapeo SYMBOL -> ENTREZ
  map_to_entrez <- function(genes_symbol) {
    out <- bitr(genes_symbol, fromType = "SYMBOL",
                toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    unique(out$ENTREZID)
  }
  
  summarize_terms <- function(enrich_obj, k = k_top_terms, term_col = "Description") {
    if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
      return(list(terms = NA_character_, padj = NA_real_))
    }
    df <- as.data.frame(enrich_obj) %>% arrange(p.adjust)
    terms <- head(df[[term_col]], k)
    padj  <- head(df$p.adjust, k)
    list(terms = paste(terms, collapse = " | "),
         padj  = paste(signif(padj, 3), collapse = " | "))
  }
  
  modules <- sort(unique(df_nodes$module))
  
  res_list <- purrr::map_df(modules, function(m){
    genes_m <- df_nodes %>% filter(module == m) %>% pull(gene) %>% unique()
    if (length(genes_m) < 3) {
      return(tibble(
        disease = disease, topN = topN, expanded = expanded,
        module = m, module_size = length(genes_m),
        GO_BP_top = NA_character_, GO_BP_padj = NA_character_,
        KEGG_top = NA_character_, KEGG_padj = NA_character_
      ))
    }
    
    entrez <- map_to_entrez(genes_m)
    
    # GO-BP
    go_bp <- tryCatch(
      enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, ont = "BP",
               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
               readable = TRUE),
      error = function(e) NULL
    )
    go_sum <- summarize_terms(go_bp, k = k_top_terms)
    
    # KEGG (organismo humano)
    kegg <- tryCatch(
      enrichKEGG(gene = entrez, organism = "hsa", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2),
      error = function(e) NULL
    )
    # Nota: enrichKEGG devuelve IDs KEGG; intentamos mapear a términos legibles si es posible
    if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
      kdf <- as.data.frame(kegg) %>% arrange(p.adjust)
      k_terms <- head(kdf$Description, k_top_terms)
      k_padj  <- head(kdf$p.adjust, k_top_terms)
      k_sum <- list(terms = paste(k_terms, collapse = " | "),
                    padj  = paste(signif(k_padj, 3), collapse = " | "))
    } else {
      k_sum <- list(terms = NA_character_, padj = NA_character_)
    }
    
    tibble(
      disease = disease, topN = topN, expanded = expanded,
      module = m, module_size = length(genes_m),
      GO_BP_top = go_sum$terms, GO_BP_padj = go_sum$padj,
      KEGG_top  = k_sum$terms, KEGG_padj  = k_sum$padj
    )
  })
  
  # ------ 5) Guardar resumen ------
  out_tsv <- glue("output/{disease}/prioritization/module_summary_top{topN}.tsv")
  readr::write_tsv(res_list %>% arrange(desc(module_size)), out_tsv)
  message(glue("Resumen por módulo (top-{topN}) guardado en: {out_tsv}"))
  invisible(res_list)
}

# Ejecutar para tus enfermedades (ajusta si lees dinámicamente)
for (d in c("CMD","CM")) {
  summarize_top_yatra_modules(disease = d, topN = 20)
}
