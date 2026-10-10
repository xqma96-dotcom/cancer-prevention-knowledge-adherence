# Diagnostic only: compare Holm correction on full-precision vs pre-rounded Wald p-values.
# Neither primary nor domain model definitions are changed here.
output_dir <- Sys.getenv("CEBP_OUTPUT_DIR", unset = file.path(getwd(), "output"))
res <- readRDS(file.path(output_dir, "02_clpm_results_N805.rds"))
models <- list(
  diet = res$fit_diet$fit, plant_fiber_wholegrain = res$fit_plantfiber_whole$fit,
  plant_fiber_fruitveg = res$fit_plantfiber_fruveg$fit,
  obesity = res$fit_obesity$fit, physical_activity = res$fit_pa$fit,
  alcohol = res$fit_alcohol$fit, smoking = res$fit_smoking$fit
)
wald_p <- function(fit, label) {
  pe <- lavaan::parameterEstimates(fit, standardized = FALSE, ci = FALSE)
  val <- pe$pvalue[pe$op == "~" & !is.na(pe$label) & pe$label == label]
  if (length(val) != 1) stop("Expected one path for label ", label)
  val
}
for (direction in c("d", "b")) {
  raw <- vapply(models, wald_p, numeric(1), label = direction)
  out <- data.frame(domain = names(raw), direction = if (direction == "d") "Knowledge-to-adherence" else "Adherence-to-knowledge",
                    raw_p = raw, legacy_p_3dp = round(raw, 3),
                    holm_raw = p.adjust(raw, "holm"),
                    holm_legacy = p.adjust(round(raw, 3), "holm"))
  print(out, row.names = FALSE)
  write.csv(out, file.path(output_dir, paste0("CHECK_Holm_precision_", direction, ".csv")), row.names = FALSE)
}
