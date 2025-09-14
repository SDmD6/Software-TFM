# =============== GENERAR INFORME (robusto a rutas) ===============
suppressPackageStartupMessages({
  library(rmarkdown); library(glue); library(fs); library(readr)
})

# --- Descubrir ruta de este script y del proyecto ---
get_script_path <- function() {
  sp <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (!is.na(sp) && nzchar(sp)) return(sp)
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    sp2 <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(sp2)) return(normalizePath(sp2))
  }
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg)) return(normalizePath(sub("^--file=", "", file_arg[1])))
  normalizePath("scripts/generate_report.R", mustWork = FALSE)
}

script_path <- get_script_path()
script_dir  <- dirname(script_path)
proj_dir    <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

# --- Parámetros por defecto (ahora que proj_dir ya existe) ---
pairs_file <- file.path(proj_dir, "data", "disease_pairs.tsv")
if (!file.exists(pairs_file)) {
  stop("No existe: ", pairs_file, "\nCrea 'data/disease_pairs.tsv' con columnas Disease y HPO_ID.")
}
pairs <- readr::read_tsv(pairs_file, show_col_types = FALSE)
if (!all(c("Disease","HPO_ID") %in% names(pairs))) {
  stop("El archivo 'disease_pairs.tsv' debe contener columnas: Disease y HPO_ID.")
}
if (nrow(pairs) < 2) stop("Se requieren al menos dos enfermedades en 'disease_pairs.tsv'.")

d1 <- pairs$Disease[1]; d2 <- pairs$Disease[2]
base_out <- "output"
base_data <- "data"
cmp_dir  <- file.path(base_out, "comparison")

# --- Candidatos para localizar la plantilla ---
candidates <- c(
  file.path(proj_dir, "report_template.Rmd"),
  file.path(getwd(),  "report_template.Rmd"),
  file.path(script_dir, "report_template.Rmd"),
  file.path(proj_dir, "reports", "report_template.Rmd"),
  file.path(proj_dir, "templates", "report_template.Rmd")
)
input_rmd <- candidates[file.exists(candidates)][1]
if (is.na(input_rmd)) {
  stop("No se encuentra 'report_template.Rmd'. Rutas probadas:\n",
       paste0(" - ", candidates, collapse = "\n"))
}

# --- Carpeta de salida + nombre del HTML ---
outdir  <- file.path(proj_dir, "reports")
fs::dir_create(outdir)
outfile <- file.path(outdir, glue("{d1}_{d2}_report.html"))

# --- Render coherente con YAML ---
# Asegúrate de que el YAML de la plantilla declare exactamente estos params:
# params:
#   project_dir: "."
#   disease_1: "CMD"
#   disease_2: "CM"
#   base_out: "output"
#   base_data: "data"
#   comparison_dir: "output/comparison"
#   use_absolute_paths: true

# --- Render coherente con YAML ---
old_root <- getOption("knitr.root.dir")
options(knitr.root.dir = proj_dir)
on.exit(options(knitr.root.dir = old_root), add = TRUE)

rmarkdown::render(
  input       = input_rmd,
  output_file = outfile,
  params = list(
    project_dir        = proj_dir,
    disease1           = d1,
    disease2           = d2,
    base_out           = base_out,
    base_data          = base_data,
    comparison_dir     = cmp_dir,
    use_absolute_paths = TRUE
  ),
  envir = new.env()
)

message(glue("Informe generado: {outfile}"))
