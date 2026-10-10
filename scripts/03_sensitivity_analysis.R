# Alternative adherence scores, retention-weighted and multigroup sensitivity analyses.


library(dplyr)
library(lavaan)
library(flextable)

OUTPUT_DIR <- Sys.getenv("CEBP_OUTPUT_DIR", unset = file.path(getwd(), "output"))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)

format_p3 <- function(x) {
  ifelse(is.na(x), NA_character_,
         ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
}

df <- readRDS(file.path(OUTPUT_DIR, "cohort_clean_N805.rds"))
cat("Loaded N =", nrow(df), "\n")

ipw_weights <- readRDS(file.path(OUTPUT_DIR, "ipw_weights_N805.rds"))
df <- df %>% left_join(ipw_weights, by = "user_id")
cat("IPW field present:", "ipw" %in% names(df),
    " | missing weights:", sum(is.na(df$ipw)), "\n")

clpm_res <- readRDS(file.path(OUTPUT_DIR, "02_clpm_results_N805.rds"))


clpm_v3_syntax <- "
  z_P_T1 ~ c*z_P_T0 + d*z_K_T0
  z_K_T1 ~ a*z_K_T0 + b*z_P_T0

  z_P_T1 + z_K_T1 ~ z_Age + Sex_Code + z_Edu + z_Income +
                     z_SES_T0 + z_Days + Work_Status + Chronic_Disease

  z_P_T0 + z_K_T0 ~ z_Age + Sex_Code + z_Edu + z_Income +
                     z_SES_T0 + Work_Status + Chronic_Disease

  z_K_T0 ~~ z_P_T0
  z_K_T1 ~~ z_P_T1
"

prep_and_fit <- function(data, p_T0, p_T1, label, use_weights = FALSE) {
  message("▶ Fitting: ", label, "  (N=", nrow(data), ")")
  d <- data %>% mutate(
    z_K_T0 = as.numeric(scale(K_Total_T0)),
    z_K_T1 = as.numeric(scale(K_Total_T1)),
    z_P_T0 = as.numeric(scale({{ p_T0 }})),
    z_P_T1 = as.numeric(scale({{ p_T1 }})),
    z_Age    = as.numeric(scale(age)),
    z_Edu    = as.numeric(scale(bs_edu)),
    z_Income = as.numeric(scale(familyIncome)),
    z_SES_T0 = as.numeric(scale(SES_3_T0)),
    z_Days   = as.numeric(scale(Days_Diff)),
    Sex_Code        = sex,
    Work_Status     = oe_m_current_T0,
    Chronic_Disease = mh_chronic_T0
  )
  if (use_weights) {
    if (!"ipw" %in% names(d) || any(is.na(d$ipw))) {
      message("  ✗ ipw column missing or contains NA, cannot fit weighted model")
      return(list(fit = NULL, label = label, data = d))
    }
  }

  fit <- tryCatch(
    if (use_weights) {
      sem(clpm_v3_syntax, data = d, estimator = "MLR", missing = "fiml",
          sampling.weights = "ipw")
    } else {
      sem(clpm_v3_syntax, data = d, estimator = "MLR", missing = "fiml")
    },
    error = function(e) { message("  ✗ ", e$message); NULL }
  )
  list(fit = fit, label = label, data = d)
}

extract_sens_row <- function(fit_obj, label, n) {
  if (is.null(fit_obj$fit)) {
    return(data.frame(Model = label, N = n,
                      KP_beta = NA, KP_CI = NA, KP_p = NA,
                      PK_beta = NA, PK_p = NA,
                      CFI = NA, RMSEA = NA, SRMR = NA, stringsAsFactors = FALSE))
  }
  fit <- fit_obj$fit
  pe <- parameterEstimates(fit, standardized = FALSE, ci = FALSE)
  ss <- standardizedSolution(fit, type = "std.all", se = TRUE, ci = TRUE, level = 0.95)
  get <- function(lbl) {
    r <- pe[!is.na(pe$label) & pe$label == lbl & pe$op == "~", ]
    if (nrow(r) == 0) return(list(beta = NA, ci = NA, p = NA))
    sr <- ss[ss$lhs == r$lhs[1] & ss$op == r$op[1] & ss$rhs == r$rhs[1], ]
    if (nrow(sr) != 1) stop("Expected exactly one standardized sensitivity path: ", lbl)
    list(beta = round(sr$est.std[1], 3),
         ci   = sprintf("[%.3f, %.3f]", sr$ci.lower[1], sr$ci.upper[1]),
         p    = round(r$pvalue[1], 3))
  }
  d_path <- get("d"); b_path <- get("b")
  fi <- fitMeasures(fit, c("cfi.robust","rmsea.robust","srmr"))
  data.frame(Model   = label, N = n,
             KP_beta = d_path$beta, KP_CI = d_path$ci, KP_p = d_path$p,
             PK_beta = b_path$beta, PK_p  = b_path$p,
             CFI     = round(fi["cfi.robust"],  3),
             RMSEA   = round(fi["rmsea.robust"], 3),
             SRMR    = round(fi["srmr"], 3), stringsAsFactors = FALSE)
}


cat("\n===== S1: WCRF version comparison =====\n")

fit_wcrf7 <- prep_and_fit(df, wcrf7_T0, wcrf7_T1, "K_Total × WCRF-7 (no smoking)")
fit_wcrf8_ref <- clpm_res$fit_main
fit_wcrf9 <- prep_and_fit(df, wcrf9_T0, wcrf9_T1, "K_Total × WCRF-9 (+hot food)")

row_wcrf7 <- extract_sens_row(fit_wcrf7, "WCRF-7 (excl. smoking, N=805)", 805)
row_wcrf8 <- extract_sens_row(fit_wcrf8_ref, "WCRF-8 (main analysis, N=805)", 805)
row_wcrf9 <- extract_sens_row(fit_wcrf9, "WCRF-9 (+hot food, N=805)", 805)

tbl_wcrf <- bind_rows(row_wcrf7, row_wcrf8, row_wcrf9)
print(tbl_wcrf[, c("Model","KP_beta","KP_CI","KP_p","CFI","RMSEA")], row.names = FALSE)


cat("\n===== S2: Follow-up ≥90 days subgroup =====\n")

df_90 <- df %>% filter(Days_Diff >= 90)
cat("N (Days_Diff ≥ 90):", nrow(df_90), "| Mean interval:", round(mean(df_90$Days_Diff), 1), "days\n")

fit_90d <- prep_and_fit(df_90, total_score1_T0, total_score1_T1, paste0("K_Total × WCRF-8, ≥90d (N=", nrow(df_90), ")"))
row_90d <- extract_sens_row(fit_90d, paste0("Interval ≥90d (N=", nrow(df_90), ")"), nrow(df_90))


cat("\n=====  Multi-group CLPM (Sex, Age, Edu, Income, Chronic) =====\n")

age_median <- median(df$age, na.rm = TRUE)

df_mgsem <- df %>%
  mutate(
    z_K_T0 = as.numeric(scale(K_Total_T0)),
    z_K_T1 = as.numeric(scale(K_Total_T1)),
    z_P_T0 = as.numeric(scale(total_score1_T0)),
    z_P_T1 = as.numeric(scale(total_score1_T1)),

    z_Age    = as.numeric(scale(age)),
    z_Edu    = as.numeric(scale(bs_edu)),
    z_Income = as.numeric(scale(familyIncome)),
    z_SES_T0 = as.numeric(scale(SES_3_T0)),
    z_Days   = as.numeric(scale(Days_Diff)),
    Sex_Code        = sex,
    Work_Status     = oe_m_current_T0,
    Chronic_Disease = mh_chronic_T0,

    Sex_Group = factor(sex, levels = c(1, 0), labels = c("Male", "Female")),
    Age_Group = factor(ifelse(age >= age_median, 1, 0), levels = c(1, 0), labels = c(sprintf("\u2265 %d Years", age_median), sprintf("< %d Years", age_median))),
    Edu_Group = factor(ifelse(bs_edu >= 2, 1, 0), levels = c(1, 0), labels = c("\u2265 Bachelor's", "< Bachelor's")),
    Inc_Group = factor(ifelse(familyIncome >= 4, 1, 0), levels = c(1, 0), labels = c("\u2265 100k CNY", "< 100k CNY")),
    Chronic_Group = factor(mh_chronic_T0, levels = c(1, 0), labels = c("Yes", "No"))
  )

run_mgsem_subgroup <- function(group_var, group_name_str, exclude_cov) {
  cat("\nRunning MG-SEM for:", group_name_str, "\n")

  all_covs <- c("z_Age", "Sex_Code", "z_Edu", "z_Income", "z_SES_T0", "z_Days", "Work_Status", "Chronic_Disease")
  all_covs <- setdiff(all_covs, exclude_cov)

  cov_T1 <- paste(all_covs, collapse = " + ")
  cov_T0 <- paste(setdiff(all_covs, "z_Days"), collapse = " + ")

  mgsem_syntax_free <- sprintf("
    z_P_T1 ~ z_P_T0 + z_K_T0
    z_K_T1 ~ z_K_T0 + z_P_T0
    z_P_T1 + z_K_T1 ~ %s
    z_P_T0 + z_K_T0 ~ %s
    z_K_T0 ~~ z_P_T0
    z_K_T1 ~~ z_P_T1
  ", cov_T1, cov_T0)

  mgsem_syntax_cons <- sprintf("
    z_P_T1 ~ z_P_T0 + c(d_eq, d_eq)*z_K_T0
    z_K_T1 ~ z_K_T0 + z_P_T0
    z_P_T1 + z_K_T1 ~ %s
    z_P_T0 + z_K_T0 ~ %s
    z_K_T0 ~~ z_P_T0
    z_K_T1 ~~ z_P_T1
  ", cov_T1, cov_T0)

  fit_free <- tryCatch(sem(mgsem_syntax_free, data = df_mgsem, group = group_var, estimator = "MLR", missing = "fiml"), error = function(e) NULL)
  fit_cons <- tryCatch(sem(mgsem_syntax_cons, data = df_mgsem, group = group_var, estimator = "MLR", missing = "fiml"), error = function(e) NULL)

  if (is.null(fit_free) || is.null(fit_cons)) return(NULL)

  lrt <- lavTestLRT(fit_free, fit_cons, method = "satorra.bentler.2010")
  delta_p <- as.numeric(lrt[2, "Pr(>Chisq)"])

  pe_raw <- parameterEstimates(
    fit_free,
    standardized = TRUE,
    ci = TRUE
  )

  pe_std <- standardizedSolution(
    fit_free,
    type = "std.all",
    se = TRUE,
    ci = TRUE
  )

  grp_labels <- lavInspect(fit_free, "group.label")

  extract_p <- function(grp_idx, lhs, rhs) {

    r <- pe_std[
      pe_std$lhs == lhs &
        pe_std$rhs == rhs &
        pe_std$op == "~" &
        pe_std$group == grp_idx,
    ]

    if (nrow(r) != 1) {
      stop("Expected exactly one matching path.")
    }

    sprintf(
      "%.3f [%.3f, %.3f]",
      r$est.std[1],
      r$ci.lower[1],
      r$ci.upper[1]
    )
  }

  extract_pval <- function(grp_idx, lhs, rhs) {

    p <- pe_raw[
      pe_raw$lhs == lhs &
        pe_raw$rhs == rhs &
        pe_raw$op == "~" &
        pe_raw$group == grp_idx,
      "pvalue"
    ]

    format_p3(p)
  }

  data.frame(
    Subgroup   = c(group_name_str, ""),
    Strata     = grp_labels,
    AR_K       = sapply(1:2, function(g) extract_p(g, "z_K_T1", "z_K_T0")),
    AR_P       = sapply(1:2, function(g) extract_p(g, "z_P_T1", "z_P_T0")),
    KP_beta_ci = sapply(1:2, function(g) extract_p(g, "z_P_T1", "z_K_T0")),
    KP_p       = sapply(1:2, function(g) extract_pval(g, "z_P_T1", "z_K_T0")),
    PK_beta_ci = sapply(1:2, function(g) extract_p(g, "z_K_T1", "z_P_T0")),
    PK_p       = sapply(1:2, function(g) extract_pval(g, "z_K_T1", "z_P_T0")),
    LRT_p      = c(format_p3(delta_p), ""),
    stringsAsFactors = FALSE
  )
}

final_mgsem_table <- bind_rows(
  run_mgsem_subgroup("Sex_Group", "Sex", "Sex_Code"),
  run_mgsem_subgroup("Age_Group", "Age", "z_Age"),
  run_mgsem_subgroup("Edu_Group", "Education", "z_Edu"),
  run_mgsem_subgroup("Inc_Group", "Household Income", "z_Income"),
  run_mgsem_subgroup("Chronic_Group", "Chronic Disease", "Chronic_Disease")
)

cat("\n===== S4: IPW-weighted primary model =====\n")

fit_ipw <- prep_and_fit(df, total_score1_T0, total_score1_T1,
                        label = "K_Total × WCRF-8 (IPW-weighted)",
                        use_weights = TRUE)

row_ipw_weighted   <- extract_sens_row(fit_ipw, "IPW-weighted (N=805)", 805)
row_ipw_unweighted <- extract_sens_row(list(fit = clpm_res$fit_main$fit),
                                       "Unweighted (main analysis, N=805)", 805)

tbl_ipw <- bind_rows(row_ipw_unweighted, row_ipw_weighted)
print(tbl_ipw[, c("Model","KP_beta","KP_CI","KP_p","PK_beta","PK_p","CFI","RMSEA")],
      row.names = FALSE)

tbl_sens_full <- bind_rows(
  data.frame(Model = "WCRF version comparison", N = NA, KP_beta=NA, KP_CI=NA, KP_p=NA, PK_beta=NA, PK_p=NA, CFI=NA, RMSEA=NA, SRMR=NA, stringsAsFactors=FALSE),
  tbl_wcrf,
  data.frame(Model = "Follow-up interval subgroup", N = NA, KP_beta=NA, KP_CI=NA, KP_p=NA, PK_beta=NA, PK_p=NA, CFI=NA, RMSEA=NA, SRMR=NA, stringsAsFactors=FALSE),
  row_90d
)

ft_S7 <- tbl_sens_full %>%
  flextable() %>%
  set_caption(caption = "Sensitivity analyses of WCRF score versions and follow-up interval.") %>%
  set_header_labels(
    Model = "Model / Subgroup", N = "N",
    KP_beta = "K\u2192P (\u03b2)", KP_CI = "K\u2192P 95% CI", KP_p = "K\u2192P p",
    PK_beta = "P\u2192K (\u03b2)", PK_p = "P\u2192K p"
  ) %>%
  bold(i = ~ !is.na(KP_p) & KP_p < 0.05, j = c("KP_beta","KP_p")) %>%
  bold(i = ~ Model %in% c("WCRF version comparison", "Follow-up interval subgroup"), bold = TRUE) %>%
  bg(i = ~ grepl("main analysis", Model), bg = "#FFF9C4", part = "body") %>%
  add_footer_lines("Models used the primary CLPM covariate specification with MLR. WCRF-7 excludes smoking; WCRF-8 is the primary adherence score; WCRF-9 adds hot-food avoidance. Coefficients and 95% CIs are standardized; path p-values are original-parameter Wald tests. Bold: p < 0.05.") %>%
  fontsize(size = 8, part = "footer") %>% italic(part = "footer") %>% autofit()

flextable::save_as_docx(ft_S7, path = file.path(OUTPUT_DIR, "WCRF_Score_and_Followup_Sensitivity.docx"))
cat("\nSaved: Sensitivity results saved\n")

print(ft_S7)

ft_S8 <- tbl_ipw %>%
  flextable() %>%
  set_caption(caption = paste0(
    "Sensitivity analysis using inverse probability ",
    "weighting (IPW) to account for differential attrition."
  )) %>%
  set_header_labels(
    Model = "Model", N = "N",
    KP_beta = "K\u2192P (\u03b2)", KP_CI = "K\u2192P 95% CI", KP_p = "K\u2192P p",
    PK_beta = "P\u2192K (\u03b2)", PK_p = "P\u2192K p"
  ) %>%
  bold(i = ~ !is.na(KP_p) & KP_p < 0.05, j = c("KP_beta","KP_p")) %>%
  add_footer_lines(paste0(
    "Note: Retention probabilities were estimated via logistic regression using ",
    "measured baseline covariates, total knowledge, and WCRF-8 adherence. Weighted and ",
    "unweighted models both used the primary covariate specification (MLR estimator; ",
    "FIML specified for missing data). Standardized coefficients and standardized 95% CIs are shown; p-values are original-parameter Wald tests. Bold indicates p < 0.05."
  )) %>%
  fontsize(size = 8, part = "footer") %>% italic(part = "footer") %>% autofit()

flextable::save_as_docx(ft_S8, path = file.path(OUTPUT_DIR, "IPW_CLPM_Sensitivity.docx"))
cat("\nSaved: IPW sensitivity results saved\n")

print(ft_S8)

ft_mgsem <- final_mgsem_table %>%
  flextable() %>%
  set_caption(caption = "Multi-group CLPM across demographic and health subgroups.") %>%
  set_header_labels(
    Subgroup   = "Subgroup", Strata = "Strata",
    AR_K       = "K \u2192 K (\u03B2 [95% CI])", AR_P = "P \u2192 P (\u03B2 [95% CI])",
    KP_beta_ci = "K \u2192 P (\u03B2 [95% CI])", KP_p = "p-value",
    PK_beta_ci = "P \u2192 K (\u03B2 [95% CI])", PK_p = "p-value",
    LRT_p      = "\u0394\u03c7\u00b2 p-value"
  ) %>%
  merge_v(j = c("Subgroup", "LRT_p")) %>%
  theme_vanilla() %>% align(align = "center", part = "all") %>% align(j = 1:2, align = "left", part = "body") %>%
  bold(i = ~ grepl("<", KP_p) | (!is.na(as.numeric(KP_p)) & as.numeric(KP_p) < 0.05), j = c("KP_beta_ci", "KP_p"), bold = TRUE) %>%
  bold(i = ~ grepl("<", PK_p) | (!is.na(as.numeric(PK_p)) & as.numeric(PK_p) < 0.05), j = c("PK_beta_ci", "PK_p"), bold = TRUE) %>%
  add_footer_lines("Note: Models estimated using MLR estimator (FIML). \u0394\u03c7\u00b2 p-value represents the Satorra-Bentler scaled LRT comparing a free model to a model where the K\u2192P path is constrained to be equal across strata. Age was dichotomized at the median.") %>%
  fontsize(size = 9, part = "footer") %>% italic(part = "footer") %>% autofit()

flextable::save_as_docx(ft_mgsem, path = file.path(OUTPUT_DIR, "Multigroup_CLPM_Subgroups.docx"))
cat("Saved: Multigroup estimates saved\n")

print(ft_mgsem)
saveRDS(
  list(
    fit_wcrf7 = fit_wcrf7, fit_wcrf9 = fit_wcrf9, tbl_wcrf = tbl_wcrf,
    fit_90d = fit_90d, row_90d = row_90d,
    mgsem_results = final_mgsem_table,
    fit_ipw = fit_ipw, tbl_ipw = tbl_ipw
  ),
  file.path(OUTPUT_DIR, "03_sensitivity_results_N805.rds")
)

cat("Saved: All objects saved → 03_sensitivity_results_N805.rds\n")
