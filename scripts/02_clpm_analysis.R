# Primary and domain-specific cross-lagged panel models.
# Standardized estimates/CIs are extracted separately from unstandardized Wald p-values.


library(dplyr)
library(lavaan)
library(flextable)
library(ggplot2)

OUTPUT_DIR <- Sys.getenv("CEBP_OUTPUT_DIR", unset = file.path(getwd(), "output"))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

df <- readRDS(file.path(OUTPUT_DIR, "cohort_clean_N805.rds"))
cat("Loaded N =", nrow(df), "\n")

set.seed(2026)


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


run_clpm <- function(data, k_var, p_var, label,
                     p_T0_vec = NULL, p_T1_vec = NULL ,
                     use_weights = FALSE) {

  message("\n▶ Running: ", label)

  k_cols <- paste0(k_var, c("_T0", "_T1"))
  missing_k <- k_cols[!k_cols %in% names(data)]
  if (length(missing_k) > 0) {
    message("  ✗ Missing K columns: ", paste(missing_k, collapse = ", "))
    return(list(fit = NULL, label = label, k_var = k_var, p_var = p_var))
  }

  if (is.null(p_T0_vec)) {
    p_cols <- paste0(p_var, c("_T0", "_T1"))
    missing_p <- p_cols[!p_cols %in% names(data)]
    if (length(missing_p) > 0) {
      message("  ✗ Missing P columns: ", paste(missing_p, collapse = ", "))
      return(list(fit = NULL, label = label, k_var = k_var, p_var = p_var))
    }
  }

  df_model <- data %>%
    mutate(
      z_K_T0 = as.numeric(scale(.data[[paste0(k_var, "_T0")]])),
      z_K_T1 = as.numeric(scale(.data[[paste0(k_var, "_T1")]])),
      z_P_T0 = if (!is.null(p_T0_vec)) as.numeric(scale(p_T0_vec))
      else as.numeric(scale(.data[[paste0(p_var, "_T0")]])),
      z_P_T1 = if (!is.null(p_T1_vec)) as.numeric(scale(p_T1_vec))
      else as.numeric(scale(.data[[paste0(p_var, "_T1")]])),
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
    if (!"ipw" %in% names(df_model) || any(is.na(df_model$ipw))) {
      message("  ✗ ipw column missing or contains NA, cannot fit weighted model")
      return(list(fit = NULL, label = label, k_var = k_var, p_var = p_var))
    }
  }

  fit <- tryCatch(
    if (use_weights) {
      sem(clpm_v3_syntax, data = df_model,
          estimator = "MLR", missing = "fiml",
          sampling.weights = "ipw")
    } else {
      sem(clpm_v3_syntax, data = df_model,
          estimator = "MLR", missing = "fiml")
    },
    error = function(e) { message("  ✗ sem() failed: ", e$message); NULL }
  )

  if (!is.null(fit) && !lavInspect(fit, "converged"))
    warning("  ⚠ Did not converge: ", label)

  list(fit = fit, label = label, k_var = k_var, p_var = p_var,weighted = use_weights)
}


extract_row <- function(fit_obj, domain, k_label, p_label) {

  empty <- data.frame(
    Domain = domain, K_Node = k_label, P_Node = p_label,
    AR_K = NA_real_, AR_P = NA_real_,
    KP_beta = NA_real_, KP_CI = NA_character_, KP_p = NA_real_,
    PK_beta = NA_real_, PK_CI = NA_character_, PK_p = NA_real_,
    CFI = NA_real_, RMSEA = NA_real_, SRMR = NA_real_,
    stringsAsFactors = FALSE
  )

  if (is.null(fit_obj) || is.null(fit_obj$fit)) return(empty)

  fit <- fit_obj$fit
  pe <- parameterEstimates(fit, standardized = FALSE, ci = FALSE)
  ss <- standardizedSolution(fit, type = "std.all", se = TRUE, ci = TRUE, level = 0.95)

  get_path <- function(lbl) {
    row <- pe[!is.na(pe$label) & pe$label == lbl & pe$op == "~", ]
    if (nrow(row) == 0) return(list(beta = NA_real_, ci = NA_character_, p = NA_real_))
    sr <- ss[ss$lhs == row$lhs[1] & ss$op == row$op[1] & ss$rhs == row$rhs[1], ]
    if (nrow(sr) != 1) stop("Expected exactly one standardized CLPM path: ", lbl)
    list(
      beta = round(sr$est.std[1], 3),
      ci   = sprintf("[%.3f, %.3f]", sr$ci.lower[1], sr$ci.upper[1]),
      p    = as.numeric(row$pvalue[1])  # Preserve precision until after multiplicity adjustment.
    )
  }

  a <- get_path("a"); c_ <- get_path("c")
  d <- get_path("d"); b  <- get_path("b")

  fi <- fitMeasures(fit, c("cfi.robust", "rmsea.robust", "srmr"))

  data.frame(
    Domain  = domain, K_Node = k_label, P_Node = p_label,
    AR_K    = a$beta,  AR_P   = c_$beta,
    KP_beta = d$beta,  KP_CI  = d$ci,   KP_p   = d$p,
    PK_beta = b$beta,  PK_CI  = b$ci,   PK_p   = b$p,
    CFI     = round(fi["cfi.robust"],  3),
    RMSEA   = round(fi["rmsea.robust"], 3),
    SRMR    = round(fi["srmr"], 3),
    stringsAsFactors = FALSE
  )
}


cat("\n===== MAIN CLPM: K_Total × WCRF-8 =====\n")

fit_main <- run_clpm(df, "K_Total", "total_score1",
                     label = "K_Total × WCRF-8 (main)")

cat("\n--- Fit indices ---\n")
print(fitMeasures(fit_main$fit,
                  c("cfi.robust","tli.robust","rmsea.robust",
                    "rmsea.ci.lower.robust","rmsea.ci.upper.robust","srmr")))

cat("\n--- Key paths (standardised β) ---\n")
parameterEstimates(fit_main$fit, standardized = TRUE) %>%
  filter(label %in% c("a","b","c","d")) %>%
  mutate(path = recode(label,
                       a = "AR_K (K_T0→K_T1)", b = "CL_PK (P_T0→K_T1)",
                       c = "AR_P (P_T0→P_T1)", d = "CL_KP (K_T0→P_T1)")) %>%
  select(path, std.all, pvalue) %>% print()


df <- df %>%
  mutate(
    P_Diet_T0 = fruvege_score_T0 + wholflour_score_T0 +
      fastfood_score_T0 + redmeat_score_T0 + soda_score_T0,
    P_Diet_T1 = fruvege_score_T1 + wholflour_score_T1 +
      fastfood_score_T1 + redmeat_score_T1 + soda_score_T1
  )

cat("\nP_Diet composite range T0:",
    min(df$P_Diet_T0, na.rm=TRUE), "–", max(df$P_Diet_T0, na.rm=TRUE),
    "| mean:", round(mean(df$P_Diet_T0, na.rm=TRUE), 2), "\n")


fit_diet <- run_clpm(df, "K_Diet", "P_Diet",
                     label = "Diet composite (K_Diet × P_Diet)",
                     p_T0_vec = df$P_Diet_T0,
                     p_T1_vec = df$P_Diet_T1)

fit_plantfiber_whole  <- run_clpm(df, "K_PlantFiber", "wholflour_score",
                                  label = "  PlantFiber × WholGrain")
fit_plantfiber_fruveg <- run_clpm(df, "K_PlantFiber", "fruvege_score",
                                  label = "  PlantFiber × FruVeg")


fit_obesity <- run_clpm(df, "K_Obesity", "BMI_score",    label = "Obesity/BMI")
fit_pa      <- run_clpm(df, "K_PA",      "pa_score",     label = "Physical Activity")
fit_alcohol <- run_clpm(df, "K_Alcohol", "alcohol_score",label = "Alcohol")
fit_smoking <- run_clpm(df, "K_Smoke",   "smoke_score",  label = "Smoking")


row_main    <- extract_row(fit_main,             "K_Total × WCRF-8 (main)",        "K_Total",      "total_score1")
row_diet    <- extract_row(fit_diet,             "Diet composite",                  "K_Diet",       "P_Diet (0-4)")
row_pf_wh   <- extract_row(fit_plantfiber_whole, "  ↳ PlantFiber × WholGrain",     "K_PlantFiber", "wholflour_score")
row_pf_fv   <- extract_row(fit_plantfiber_fruveg,"  ↳ PlantFiber × FruVeg",        "K_PlantFiber", "fruvege_score")
row_obesity <- extract_row(fit_obesity,          "Obesity / BMI",                   "K_Obesity",    "BMI_score")
row_pa      <- extract_row(fit_pa,               "Physical Activity",               "K_PA",         "pa_score")
row_alcohol <- extract_row(fit_alcohol,          "Alcohol",                         "K_Alcohol",    "alcohol_score")
row_smoking <- extract_row(fit_smoking,          "Smoking",                         "K_Smoke",      "smoke_score")

table3_data <- bind_rows(
  row_main,
  row_diet,
  row_pf_wh,
  row_pf_fv,
  row_obesity,
  row_pa,
  row_alcohol,
  row_smoking
)

cat("\n===== CLPM estimate preview =====\n")
print(table3_data[, c("Domain","AR_K","AR_P","KP_beta","KP_CI","KP_p","PK_beta","PK_p")],
      row.names = FALSE)


is_main <- table3_data$Domain == "K_Total × WCRF-8 (main)"

table3_data <- table3_data %>%
  mutate(
    KP_p_holm = NA_real_,
    KP_q_fdr  = NA_real_,
    PK_p_holm = NA_real_,
    PK_q_fdr  = NA_real_
  )

table3_data$KP_p_holm[!is_main] <- p.adjust(table3_data$KP_p[!is_main], method = "holm")
table3_data$KP_q_fdr[!is_main]  <- p.adjust(table3_data$KP_p[!is_main], method = "BH")

table3_data$PK_p_holm[!is_main] <- p.adjust(table3_data$PK_p[!is_main], method = "holm")
table3_data$PK_q_fdr[!is_main]  <- p.adjust(table3_data$PK_p[!is_main], method = "BH")

cat("\n===== Multiplicity adjustments: seven domain models K-to-P (primary excluded)=====\n")
print(table3_data[, c("Domain","KP_beta","KP_p","KP_p_holm","KP_q_fdr")], row.names = FALSE)

cat("\n===== Multiplicity adjustments: seven domain models P-to-K (primary excluded)=====\n")
print(table3_data[, c("Domain","PK_beta","PK_p","PK_p_holm","PK_q_fdr")], row.names = FALSE)

sig_index <- which(
  (is_main & !is.na(table3_data$KP_p) & table3_data$KP_p < 0.05) |
    (!is_main & !is.na(table3_data$KP_p_holm) & table3_data$KP_p_holm < 0.05)
)

sig_index_pk <- which(
  (is_main & !is.na(table3_data$PK_p) & table3_data$PK_p < 0.05) |
    (!is_main & !is.na(table3_data$PK_p_holm) & table3_data$PK_p_holm < 0.05)
)

tbl3_display <- table3_data %>%
  mutate(
    KP_p_display = ifelse(
      is.na(KP_p_holm),
      sprintf("%.3f", KP_p),
      sprintf("%.3f (%.3f)", KP_p, KP_p_holm)
    ),
    PK_p_display = ifelse(
      is.na(PK_p_holm),
      sprintf("%.3f", PK_p),
      sprintf("%.3f (%.3f)", PK_p, PK_p_holm)
    )
  ) %>%
  select(
    Domain,
    `Knowledge\nConstruct` = K_Node,
    `Behavior\nConstruct`  = P_Node,
    `AR K→K\n(β)`          = AR_K,
    `AR P→P\n(β)`          = AR_P,
    `K→P\n(β)`             = KP_beta,
    `K→P\n95% CI`          = KP_CI,
    `K→P\np (Holm)`        = KP_p_display,
    `P→K\n(β)`             = PK_beta,
    `P→K\n95% CI`          = PK_CI,
    `P→K\np (Holm)`        = PK_p_display
  )

ft3 <- flextable(tbl3_display) %>%
  set_caption(
    caption = as_paragraph(
      "Cross-lagged panel model estimates for overall and domain-specific knowledge–behavior associations."
    )
  ) %>%
  bg(j = c("K→P\n(β)", "K→P\n95% CI", "K→P\np (Holm)"),
     bg = "#DDEEFF", part = "body") %>%
  bg(j = c("K→P\n(β)", "K→P\n95% CI", "K→P\np (Holm)"),
     bg = "#BBCCEE", part = "header") %>%
  bold(i = sig_index,
       j = c("K→P\n(β)", "K→P\np (Holm)")) %>%
  bold(
    i = sig_index_pk,
    j = c("P→K\n(β)", "P→K\np (Holm)")
  ) %>%
  italic(
    i = grep("↳", table3_data$Domain),
    j = "Domain"
  ) %>%
  border(i = 1, border.bottom = officer::fp_border(color = "black", width = 1.5)) %>%
  align(align = "center", part = "header") %>%
  align(j = 4:11, align = "center", part = "body") %>%
  align(j = 1:3,  align = "left",   part = "body") %>%
  add_footer_lines(values = paste0(
    "All models used the prespecified covariate specification (MLR estimator, FIML for missing data), ",
    "controlling for age, sex, education, family income, subjective SES, follow-up interval, ",
    "employment status, and chronic disease history. ",
    "Standardized \u03b2 coefficients are reported. ",
    "AR = autoregressive path; K\u2192P = prospective association of baseline knowledge with follow-up adherence; ",
    "P\u2192K = prospective association of baseline adherence with follow-up knowledge. ",
    "Standardized 95% CIs are from lavaan::standardizedSolution; path p-values are unstandardized-parameter Wald tests. ",
    "Diet composite (P_Diet) = sum of five WCRF/AICR diet sub-scores: ",
    "fruit & vegetables, whole grains, fast food, red/processed meat, sugar-sweetened beverages (range 0\u20134). ",
    "Indented rows show knowledge sub-domain \u00d7 individual behavior pairings. ",
    "Model fit indices are reported in the model fit output. ",
    "K_Total (main model) was pre-specified as the primary hypothesis and not subjected to multiplicity correction, ",
    "as the 10 knowledge sub-domains derive from theoretical constructs of the CCAC scale rather than a data-driven partition. ",
    "The 7 domain-specific K\u2192P associations and the 7 domain-specific P\u2192K associations were each Holm-corrected as separate, ",
    "direction-specific families (family-wise); bold indicates p < .05 after Holm correction for sub-domain rows, ",
    "or uncorrected p < .05 for the K_Total row. Numbers in parentheses are Holm-corrected p-values; unparenthesized numbers are uncorrected. ",
    "See the multiplicity-adjusted estimates output for FDR-BH q-values (sensitivity analysis, both directions)."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  italic(part = "footer") %>%
  autofit()

ft3

flextable::save_as_docx(
  ft3,
  path = file.path(OUTPUT_DIR, "CLPM_Standardized_Path_Estimates.docx")
)
cat("\nSaved: the CLPM estimates saved\n")


tableS5_data <- table3_data %>%
  filter(!is_main) %>%
  select(Domain, KP_beta, KP_p, KP_p_holm, KP_q_fdr,
         PK_beta, PK_p, PK_p_holm, PK_q_fdr)

ftS5 <- flextable(tableS5_data) %>%
  colformat_double(
    j = c("KP_beta", "KP_p", "KP_p_holm", "KP_q_fdr",
          "PK_beta", "PK_p", "PK_p_holm", "PK_q_fdr"),
    digits = 3
  ) %>%
  set_header_labels(
    Domain    = "Domain",
    KP_beta   = "K\u2192P (\u03b2)",   KP_p = "K\u2192P Uncorrected p",
    KP_p_holm = "K\u2192P Holm p",     KP_q_fdr = "K\u2192P FDR-BH q",
    PK_beta   = "P\u2192K (\u03b2)",   PK_p = "P\u2192K Uncorrected p",
    PK_p_holm = "P\u2192K Holm p",     PK_q_fdr = "P\u2192K FDR-BH q"
  ) %>%
  set_caption("Multiplicity-corrected knowledge-to-adherence and adherence-to-knowledge estimates.") %>%
  bg(part = "header", bg = "#1a4a6b") %>%
  color(part = "header", color = "white") %>%
  bold(part = "header") %>%
  bold(i = which(tableS5_data$KP_q_fdr < 0.05), j = c("KP_q_fdr")) %>%
  bold(i = which(tableS5_data$PK_q_fdr < 0.05), j = c("PK_q_fdr")) %>%
  align(align = "center", part = "all") %>%
  align(j = "Domain", align = "left", part = "body") %>%
  fontsize(size = 9, part = "all") %>%
  add_footer_lines(values = paste0(
    "Note. Holm p = Holm-corrected p-value within each direction-specific family of 7 domain-specific associations ",
    "(K\u2192P and P\u2192K are corrected separately; main analysis conclusions reported in the CLPM estimates). ",
    "FDR-BH q = Benjamini-Hochberg adjusted q-value within the same respective family (sensitivity analysis). ",
    "K_Total (primary hypothesis) is excluded from both correction families and reported uncorrected in the CLPM estimates."
  )) %>%
  fontsize(size = 7, part = "footer") %>%
  italic(part = "footer") %>%
  autofit()

ftS5

save_as_docx(ftS5, path = file.path(OUTPUT_DIR, "CLPM_Multiplicity_Adjustments.docx"))
cat("Saved: CLPM_Multiplicity_Adjustments.docx saved\n")


fit_objects <- list(
  list(obj = fit_main,            label = "Main: K_Total × WCRF-8"),
  list(obj = fit_diet,            label = "Diet composite"),
  list(obj = fit_plantfiber_whole,label = "PlantFiber × WholGrain"),
  list(obj = fit_plantfiber_fruveg,label= "PlantFiber × FruVeg"),
  list(obj = fit_obesity,         label = "Obesity/BMI"),
  list(obj = fit_pa,              label = "Physical Activity"),
  list(obj = fit_alcohol,         label = "Alcohol"),
  list(obj = fit_smoking,         label = "Smoking")
)

tbl_S6 <- bind_rows(lapply(fit_objects, function(x) {
  if (is.null(x$obj$fit)) {
    return(data.frame(Model = x$label, Chi2 = NA, df = NA, p_chi2 = NA,
                      CFI = NA, RMSEA = NA, RMSEA_lo = NA, RMSEA_hi = NA,
                      SRMR = NA, stringsAsFactors = FALSE))
  }
  fi <- fitMeasures(x$obj$fit, c(
    "chisq.scaled","df.scaled","pvalue.scaled",
    "cfi.robust","rmsea.robust",
    "rmsea.ci.lower.robust","rmsea.ci.upper.robust","srmr"
  ))
  data.frame(
    Model    = x$label,
    Chi2     = round(fi["chisq.scaled"], 2),
    df       = fi["df.scaled"],
    p_chi2   = round(fi["pvalue.scaled"], 3),
    CFI      = round(fi["cfi.robust"], 3),
    RMSEA    = round(fi["rmsea.robust"], 3),
    RMSEA_lo = round(fi["rmsea.ci.lower.robust"], 3),
    RMSEA_hi = round(fi["rmsea.ci.upper.robust"], 3),
    SRMR     = round(fi["srmr"], 3),
    stringsAsFactors = FALSE
  )
}))

cat("\n=== the model fit output: Model fit indices ===\n")
print(tbl_S6, row.names = FALSE)

tbl_S6 %>%
  flextable() %>%
  set_caption(caption = "Fit indices for all CLPM models.") %>%
  set_header_labels(
    Model    = "Model",
    Chi2     = "\u03c7\u00b2 (scaled)",
    df       = "df",
    p_chi2   = "p",
    CFI      = "CFI (robust)",
    RMSEA    = "RMSEA (robust)",
    RMSEA_lo = "90% CI lower",
    RMSEA_hi = "90% CI upper",
    SRMR     = "SRMR"
  ) %>%
  autofit() %>%
  flextable::save_as_docx(
    path = file.path(OUTPUT_DIR, "CLPM_Model_Fit.docx")
  )
cat("Saved: Model fit table saved\n")


  pe_main_raw <- parameterEstimates(fit_main$fit, standardized = FALSE, ci = FALSE)
pe_main_std <- standardizedSolution(fit_main$fit, type = "std.all", se = TRUE, ci = TRUE, level = 0.95)
pe_main <- merge(
  pe_main_raw[, c("lhs", "op", "rhs", "label", "pvalue")],
  pe_main_std[, c("lhs", "op", "rhs", "est.std", "ci.lower", "ci.upper")],
  by = c("lhs", "op", "rhs"), sort = FALSE
)
pe_main$std.all <- pe_main$est.std

star_p <- function(p) {
  as.character(cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                   labels = c("***", "**", "*", "ns")))
}

paths <- pe_main %>%
  filter(op == "~", !is.na(label), label %in% c("a", "b", "c", "d")) %>%
  transmute(
    path = case_when(
      label == "a" ~ "K_T0 -> K_T1 (AR K->K)",
      label == "c" ~ "P_T0 -> P_T1 (AR P->P)",
      label == "d" ~ "K_T0 -> P_T1 (cross-lagged K->P)",
      label == "b" ~ "P_T0 -> K_T1 (cross-lagged P->K)"
    ),
    beta = round(std.all, 3),
    ci   = sprintf("[%.3f, %.3f]", ci.lower, ci.upper),
    p    = round(pvalue, 3),
    sig  = star_p(pvalue)
  )

r_T0 <- pe_main %>%
  filter(op == "~~",
         (lhs == "z_K_T0" & rhs == "z_P_T0") | (lhs == "z_P_T0" & rhs == "z_K_T0")) %>%
  transmute(path = "K_T0 <-> P_T0 (baseline r)",
            beta = round(std.all, 3),
            ci   = sprintf("[%.3f, %.3f]", ci.lower, ci.upper),
            p    = round(pvalue, 3),
            sig  = star_p(pvalue))

r_T1 <- pe_main %>%
  filter(op == "~~",
         (lhs == "z_K_T1" & rhs == "z_P_T1") | (lhs == "z_P_T1" & rhs == "z_K_T1")) %>%
  transmute(path = "K_T1 <-> P_T1 (residual r)",
            beta = round(std.all, 3),
            ci   = sprintf("[%.3f, %.3f]", ci.lower, ci.upper),
            p    = round(pvalue, 3),
            sig  = star_p(pvalue))

fig2_panelA_data <- bind_rows(paths, r_T0, r_T1)

cat("\n===== Primary CLPM figure: path estimates and significance=====\n")
print(fig2_panelA_data, row.names = FALSE)


kp_ci_num <- table3_data$KP_CI %>%
  gsub("\\[|\\]", "", .) %>%
  strsplit(", ") %>%
  lapply(as.numeric)

fig2b_data <- table3_data %>%
  mutate(
    KP_ci_lo = sapply(kp_ci_num, `[`, 1),
    KP_ci_hi = sapply(kp_ci_num, `[`, 2),
    is_main  = Domain == "K_Total × WCRF-8 (main)",
    sig      = ifelse(is_main, KP_p < 0.05, KP_p_holm < 0.05),
    label = case_when(
      Domain == "K_Total × WCRF-8 (main)"     ~ "Knowledge Total × WCRF-8 (Primary)",
      Domain == "Diet composite"               ~ "Diet composite",
      grepl("WholGrain", Domain) ~
        "    → PlantFiber × Dietary Fiber",
      grepl("FruVeg", Domain)                  ~ "    \u2192 PlantFiber \u00d7 FruVeg",
      Domain == "Obesity / BMI"                ~ "Obesity / BMI",
      Domain == "Physical Activity"            ~ "Physical Activity",
      Domain == "Alcohol"                      ~ "Alcohol",
      Domain == "Smoking"                      ~ "Smoking",
      TRUE                                     ~ Domain
    ),
    grp = case_when(
      is_main & sig  ~ "Primary, significant",
      is_main & !sig ~ "Primary, non-significant",
      !is_main & sig  ~ "Sub-domain, significant",
      !is_main & !sig ~ "Sub-domain, non-significant"
    )
  )

order_domains <- c(
  "K_Total × WCRF-8 (main)",
  "Diet composite",
  grep("WholGrain", table3_data$Domain, value = TRUE),
  grep("FruVeg",    table3_data$Domain, value = TRUE),
  table3_data %>%
    filter(!Domain %in% c("K_Total × WCRF-8 (main)", "Diet composite") &
             !grepl("WholGrain|FruVeg", Domain)) %>%
    arrange(desc(AR_P)) %>%
    pull(Domain)
)

fig2b_data <- fig2b_data %>%
  mutate(Domain = factor(Domain, levels = rev(order_domains))) %>%
  arrange(Domain)

fig2b_data <- fig2b_data %>%
  mutate(star_x = KP_ci_hi + 0.012)

fig2_panelB <- ggplot(fig2b_data, aes(x = KP_beta, y = Domain, color = grp)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = KP_ci_lo, xmax = KP_ci_hi), height = 0.18, linewidth = 0.6) +
  geom_point(aes(shape = grp), size = 3) +
  geom_text(
    data = filter(fig2b_data, sig),
    aes(x = star_x, label = "*"),
    color = "black", size = 6, vjust = 0.75, show.legend = FALSE
  ) +
  scale_y_discrete(labels = setNames(fig2b_data$label, fig2b_data$Domain)) +
  scale_color_manual(values = c(
    "Primary, significant"        = "#C2185B",
    "Primary, non-significant"    = "#C2185B",
    "Sub-domain, significant"     = "#C2185B",
    "Sub-domain, non-significant" = "#ABABAB"
  )) +
  scale_shape_manual(values = c(
    "Primary, significant"        = 18,
    "Primary, non-significant"    = 18,
    "Sub-domain, significant"     = 16,
    "Sub-domain, non-significant" = 16
  )) +
  labs(
    x = expression(paste("Standardized ", beta, " (K", symbol("\256"), "P)")),
    y = NULL,
    title = "B",
    color = NULL, shape = NULL
  ) +
  theme_classic(base_size = 11,base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 13, hjust = -0.15),
    axis.text.y = element_text(color = "black")
  )

ggsave(
  filename = file.path(OUTPUT_DIR, "CLPM_Knowledge_to_Adherence_Forest.png"),
  plot = fig2_panelB, width = 6, height = 4.2, dpi = 300, bg = "white",
  device = ragg::agg_png
)
cat("Saved: CLPM forest plot (forest plot) saved\n")
print(fig2_panelB)


saveRDS(
  list(
    fit_main             = fit_main,
    fit_diet             = fit_diet,
    fit_plantfiber_whole = fit_plantfiber_whole,
    fit_plantfiber_fruveg= fit_plantfiber_fruveg,
    fit_obesity          = fit_obesity,
    fit_pa               = fit_pa,
    fit_alcohol          = fit_alcohol,
    fit_smoking          = fit_smoking,
    table3_data          = table3_data,
    fit_index_table      = tbl_S6
  ),
  file.path(OUTPUT_DIR, "02_clpm_results_N805.rds")
)
cat("\nSaved: All objects saved → 02_clpm_results_N805.rds\n")
cat("=== 02_CLPM_analysis.R complete ===\n")
