# Compare key manuscript paths with standardized coefficients and matched confidence intervals.
# Note: displayed p-values in the manuscript use unstandardized Wald tests.
library(lavaan)
out_dir <- Sys.getenv("CEBP_OUTPUT_DIR", unset = file.path(getwd(), "output"))
cl <- readRDS(file.path(out_dir, "02_clpm_results_N805.rds"))
stopifnot(!is.null(cl$fit_main$fit), !is.null(cl$fit_diet$fit), !is.null(cl$fit_plantfiber_fruveg$fit))
extract <- function(x, direction) {
  ss <- lavaan::standardizedSolution(
    x$fit, type = "std.all", se = TRUE, ci = TRUE
  )
  ss <- ss[
    ss$op == "~" &
      !is.na(ss$label) &
      ss$label == direction, 
  ]
  
  pe <- lavaan::parameterEstimates(
    x$fit, standardized = FALSE, ci = FALSE
  )
  pe <- pe[
    pe$op == "~" &
      !is.na(pe$label) &
      pe$label == direction,
  ]
  
  if (nrow(ss) != 1L || nrow(pe) != 1L) {
    stop("Path not uniquely identified")
  }
  
  c(
    beta = ss$est.std,
    lo = ss$ci.lower,
    hi = ss$ci.upper,
    p = pe$pvalue
  )
}

expected <- list(
  main_K_to_P = list(model=cl$fit_main, path="d", beta=-0.080),
  main_P_to_K = list(model=cl$fit_main, path="b", beta=-0.003),
  diet_K_to_P = list(model=cl$fit_diet, path="d", beta=-0.102),
  fruitveg_P_to_plantfiber_K = list(model=cl$fit_plantfiber_fruveg, path="b", beta=0.084)
)
rows <- lapply(names(expected), function(nm) {
  z <- expected[[nm]]
  s <- extract(z$model,z$path)
  data.frame(path=nm, reported_beta=z$beta, extracted_beta=s[["beta"]],
             standardized_CI_lower=s[["lo"]], standardized_CI_upper=s[["hi"]],
             unadjusted_p=s[["p"]],
             agrees_at_three_decimals=(round(s[["beta"]],3) == z$beta))
})
ans <- do.call(rbind, rows)
print(ans, row.names=FALSE)
write.csv(ans, file.path(out_dir,"CHECK_standardized_coefficients.csv"),row.names=FALSE)
cat("NOTE: Domain adjusted p-values come from seven-test Holm correction and are not recalculated here.\n")
cat("NOTE: Inspect manuscript 95% CIs against standardizedSolution() output.\n")
