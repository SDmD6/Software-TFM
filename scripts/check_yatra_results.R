#Comprobaciones rápidas
library(tidyverse); library(glue)

for (d in c("CMD","CM")) {
  comb <- read_tsv(glue("output/{d}/prioritization/combined_ranking.tsv"), show_col_types = FALSE)
  stopifnot(all(c("gene","score_yatra","score_topo","score_combined") %in% names(comb)))
  cat(d, ": filas =", nrow(comb), 
      "| genes únicos =", length(unique(comb$gene)),
      "| top-1 =", comb$gene[1], "\n")
}

#Figura: correlación Yatra vs. topología 
plots_dir <- "output/comparison/plots"
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

plot_corr <- function(d) {
  comb <- read_tsv(glue("output/{d}/prioritization/combined_ranking.tsv"), show_col_types = FALSE) %>%
    drop_na(score_yatra, score_topo)
  r <- suppressWarnings(cor(comb$score_yatra, comb$score_topo, method = "spearman"))
  p <- ggplot(comb, aes(score_topo, score_yatra)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = glue("{d} · Yatra vs topología (ρ={round(r,3)})"),
         x = "Score topológico (z-degree + z-betweenness)", y = "Score Yatra")
  ggsave(glue("{plots_dir}/{d}_corr_yatra_topo.png"), p, width = 6, height = 4, dpi = 160)
}

lapply(c("CMD","CM"), plot_corr)
cat("Figuras guardadas en ", plots_dir, "\n", sep = "")

#Métrica: solapamiento top-N 
topN <- 50

get_top <- function(d, col) {
  read_tsv(glue("output/{d}/prioritization/combined_ranking.tsv"), show_col_types = FALSE) %>%
    arrange(desc(.data[[col]])) %>%
    slice_head(n = topN) %>% pull(gene)
}

top_cmd_y <- get_top("CMD","score_yatra"); top_cm_y <- get_top("CM","score_yatra")
top_cmd_t <- get_top("CMD","score_topo");  top_cm_t <- get_top("CM","score_topo")

jacc <- function(a,b) length(intersect(a,b))/length(unique(c(a,b)))

cat("Jaccard top", topN, "Yatra (CMD vs CM):", round(jacc(top_cmd_y, top_cm_y),3), "\n")
cat("Jaccard top", topN, "Topológico (CMD vs CM):", round(jacc(top_cmd_t, top_cm_t),3), "\n")

