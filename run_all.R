# Run from the root of this folder: Rscript run_all.R
# Required private inputs: data/base.sas7bdat and data/follow.sas7bdat.
# To override paths, set CEBP_DATA_DIR / CEBP_OUTPUT_DIR environment variables.
root <- normalizePath(getwd(), mustWork = TRUE)
data_dir <- Sys.getenv("CEBP_DATA_DIR", unset = file.path(root, "data"))
output_dir <- Sys.getenv("CEBP_OUTPUT_DIR", unset = file.path(root, "output"))
Sys.setenv(CEBP_DATA_DIR = data_dir, CEBP_OUTPUT_DIR = output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
required_files <- file.path(data_dir, c("base.sas7bdat", "follow.sas7bdat"))
if (!all(file.exists(required_files))) {
  stop("Required restricted source data not found: ", paste(basename(required_files[!file.exists(required_files)]), collapse = ", "),
       ". Do not upload identifiable participant data to a public repository.")
}
scripts <- c("00_data_cleaning.R", "01_descriptive_analysis.R", "02_clpm_analysis.R",
             "03_sensitivity_analysis.R", "04_network_analysis.R")
for (script in scripts) {
  message("\n===== Running ", script, " =====")
  # Global session is intentional: 01 consumes the in-memory raw/clean objects from 00.
  source(file.path(root, "scripts", script), local = FALSE, encoding = "UTF-8")
}
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
message("Completed all five stages. Compare numeric outputs with the manuscript before public release.")
