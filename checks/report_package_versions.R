# Record the versions of analysis dependencies in the local R installation.
# This does not load restricted data or execute any statistical models.
packages <- c(
  "dplyr", "haven", "lubridate", "gtsummary", "flextable", "officer",
  "ggplot2", "lavaan", "ragg", "moments", "survey", "cobalt",
  "boot", "bootnet", "qgraph", "networktools",
  "NetworkComparisonTest", "igraph", "tidyr", "patchwork",
  "gridExtra", "scales"
)
installed <- utils::installed.packages()
found <- packages %in% rownames(installed)
versions <- rep(NA_character_, length(packages))
versions[found] <- installed[packages[found], "Version"]
report <- data.frame(package = packages, version = versions,
                     installed = found, stringsAsFactors = FALSE)
output_dir <- Sys.getenv("CEBP_OUTPUT_DIR", unset = file.path(getwd(), "output"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(report, file.path(output_dir, "package_versions.csv"), row.names = FALSE)
print(report, row.names = FALSE)
