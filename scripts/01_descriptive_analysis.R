# Descriptive statistics, attrition comparisons, ceiling effects and IPW.
# Run after 00 in the same R session.


library(dplyr)
library(haven)
library(gtsummary)
library(flextable)

DATA_DIR <- Sys.getenv("CEBP_DATA_DIR", unset = file.path(getwd(), "data"))
OUTPUT_DIR <- Sys.getenv("CEBP_OUTPUT_DIR", unset = file.path(getwd(), "output"))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

required_session_objects <- c("base_raw", "cohort_clean_final", "score_knowledge", "n_base_valid")
missing_session_objects <- required_session_objects[!vapply(required_session_objects, exists, logical(1), inherits = TRUE)]
if (length(missing_session_objects)) stop("Run via run_all.R: 01 requires objects from 00: ", paste(missing_session_objects, collapse = ", "))
df <- readRDS(file.path(OUTPUT_DIR, "cohort_clean_N805.rds"))
cat("Loaded N =", nrow(df), "\n")

set.seed(2026)


df_tbl1 <- df %>%
  mutate(
    Sex          = factor(sex, levels = c(1, 0), labels = c("Male", "Female")),
    Age          = age,
    Education    = factor(bs_edu, levels = 0:3, labels = c("Primary or below", "Secondary", "Bachelor's", "Postgraduate")),
    Income       = factor(familyIncome, levels = 1:6, labels = c("<30k", "30-50k", "50-100k", "100-200k", "200-500k", ">500k")),
    SES_Score    = SES_3_T0,
    Employed     = factor(oe_m_current_T0, levels = c(1, 0), labels = c("Yes", "No")),
    Chronic_Dis  = factor(mh_chronic_T0, levels = c(0, 1), labels = c("No", "Yes")),
    Followup_day = Days_Diff
  )

set.seed(2026)
tbl1 <- df_tbl1 %>%
  select(
    Sex, Age, Education, `Annual income (CNY)` = Income,
    `Subjective SES (1-10)` = SES_Score,
    `Currently employed` = Employed,
    `Chronic disease` = Chronic_Dis,
    `Follow-up interval (days)` = Followup_day,
    `Knowledge Total Score (max:59)` = K_Total_T0,
    `Knowledge: Smoking` = K_Smoke_T0,
    `Knowledge: Alcohol` = K_Alcohol_T0,
    `Knowledge: Obesity` = K_Obesity_T0,
    `Knowledge: Physical Activity` = K_PA_T0,
    `Knowledge: Diet (composite)` = K_Diet_T0,
    `Knowledge: Occupational Exposure` = K_Occu_T0,
    `Knowledge: Infection` = K_Infection_T0,
    `Knowledge: Sun Exposure` = K_Sun_T0,
    `Knowledge: Air Pollution` = K_Air_T0,
    `Knowledge: Breastfeeding` = K_Breast_T0,
    `2018 WCRF/AICR behavior score (0-8)` = total_score1_T0
  ) %>%
  tbl_summary(
    by = Sex,
    type = list(
      starts_with("Knowledge") ~ "continuous",
      starts_with("2018 WCRF") ~ "continuous",
      `Subjective SES (1-10)` ~ "continuous"
    ),
    statistic = list(
      all_continuous() ~ "{median} [{p25}, {p75}]",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    missing = "no"
  ) %>%
  add_overall() %>%
  add_p(
    test = list(all_continuous() ~ "wilcox.test", all_categorical() ~ "chisq.test"),
    test.args = all_categorical() ~ list(simulate.p.value = TRUE),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  modify_header(label = "**Variable**") %>%
  modify_caption("**Baseline characteristics of the analytic cohort (N=805)**") %>%
  bold_labels()

footnote_text <- "Data are presented as median [IQR] for continuous variables and n (%) for categorical variables. P-values are from the Mann-Whitney U test for continuous variables, and the chi-squared test (or Fisher's exact test with Monte Carlo simulation for sparse data) for categorical variables. Non-parametric tests were applied to all continuous variables, as several knowledge subscales showed pronounced ceiling effects. Knowledge domain scores are derived from 34 items with varying item weights (theoretical score range: 0-59). WCRF = World Cancer Research Fund; AICR = American Institute for Cancer Research; SES = socioeconomic status. "
footnote_text <- paste0(
  footnote_text,
  " For knowledge subscales with a narrow possible score range (e.g., \u22643 distinct values), ",
  "median [IQR] may not fully capture distributional differences detected by the rank-based Wilcoxon test; "
)

tbl1 %>%
  as_flex_table() %>%
  delete_part(part = "footer") %>%
  add_footer_lines(values = footnote_text) %>%
  fontsize(size = 9, part = "footer") %>%
  italic(part = "footer") %>%
  flextable::save_as_docx(path = file.path(OUTPUT_DIR, "Baseline_Characteristics_by_Sex.docx"))

print(tbl1)


make_change_row <- function(df, var_T0, var_T1, label, is_total = FALSE) {
  valid_idx <- complete.cases(df[[var_T0]], df[[var_T1]])
  x0 <- df[[var_T0]][valid_idx]
  x1 <- df[[var_T1]][valid_idx]
  if (length(x0) == 0) return(NULL)

  wt <- wilcox.test(x0, x1, paired = TRUE)
  diffs <- x1 - x0
  n_total <- length(diffs)

  imp <- sprintf("%.1f%%", sum(diffs > 0) / n_total * 100)
  unc <- sprintf("%.1f%%", sum(diffs == 0) / n_total * 100)
  dec <- sprintf("%.1f%%", sum(diffs < 0) / n_total * 100)

  if (is_total) {
    t0_out <- sprintf("%.1f [%.1f, %.1f]", median(x0), quantile(x0, 0.25), quantile(x0, 0.75))
    t1_out <- sprintf("%.1f [%.1f, %.1f]", median(x1), quantile(x1, 0.25), quantile(x1, 0.75))
  } else {
    t0_out <- "-"
    t1_out <- "-"
  }

  data.frame(
    Variable      = label,
    T0_Med_IQR    = t0_out,
    T1_Med_IQR    = t1_out,
    Improved      = imp,
    Unchanged     = unc,
    Declined      = dec,
    Wilcoxon_P    = style_pvalue(wt$p.value, digits = 3),
    stringsAsFactors = FALSE
  )
}

var_list <- list(
  list("K_Total_T0", "K_Total_T1", "Knowledge total score (0-59)", TRUE),
  list("K_Smoke_T0", "K_Smoke_T1", "  Smoking", FALSE),
  list("K_Alcohol_T0", "K_Alcohol_T1", "  Alcohol", FALSE),
  list("K_Obesity_T0", "K_Obesity_T1", "  Obesity", FALSE),
  list("K_PA_T0", "K_PA_T1", "  Physical activity", FALSE),
  list("K_Diet_T0", "K_Diet_T1", "  Diet", FALSE),
  list("K_Occu_T0", "K_Occu_T1", "  Occupation", FALSE),
  list("K_Infection_T0", "K_Infection_T1", "  Infection", FALSE),
  list("K_Sun_T0", "K_Sun_T1", "  Sun exposure", FALSE),
  list("K_Air_T0", "K_Air_T1", "  Air pollution", FALSE),
  list("K_Breast_T0", "K_Breast_T1", "  Breastfeeding", FALSE),
  list("total_score1_T0", "total_score1_T1", "WCRF behavior score (0-8)", TRUE),
  list("smoke_score_T0", "smoke_score_T1", "  Smoking", FALSE),
  list("BMI_score_T0", "BMI_score_T1", "  BMI", FALSE),
  list("pa_score_T0", "pa_score_T1", "  Physical activity", FALSE),
  list("fruvege_score_T0", "fruvege_score_T1", "  Fruit & vegetables", FALSE),
  list("wholflour_score_T0", "wholflour_score_T1", "  Whole grains", FALSE),
  list("fastfood_score_T0", "fastfood_score_T1", "  Fast food", FALSE),
  list("redmeat_score_T0", "redmeat_score_T1", "  Red & processed meat", FALSE),
  list("soda_score_T0", "soda_score_T1", "  Sugar-sweetened beverage", FALSE),
  list("alcohol_score_T0", "alcohol_score_T1", "  Alcohol", FALSE)
)

rows <- lapply(var_list, function(v) {
  if (all(c(v[[1]], v[[2]]) %in% names(df))) {
    make_change_row(df, v[[1]], v[[2]], v[[3]], v[[4]])
  }
})
tbl2 <- do.call(rbind, Filter(Negate(is.null), rows))

empty_row <- function(domain_name) {
  data.frame(Variable = domain_name, T0_Med_IQR = NA, T1_Med_IQR = NA,
             Improved = NA, Unchanged = NA, Declined = NA, Wilcoxon_P = NA)
}
tbl2 <- bind_rows(empty_row("Knowledge"), tbl2[1:11, ], empty_row("Behavior"), tbl2[12:nrow(tbl2), ])

tbl2 %>%
  flextable() %>%
  set_caption(caption = "Changes in cancer prevention knowledge and behaviors from baseline to follow-up.") %>%
  set_header_labels(
    Variable   = "Variable",
    T0_Med_IQR = "Baseline Median [IQR]",
    T1_Med_IQR = "Follow-up Median [IQR]",
    Improved   = "Improved",
    Unchanged  = "Unchanged",
    Declined   = "Declined",
    Wilcoxon_P = "Wilcoxon p"
  ) %>%
  bold(i = ~ Variable %in% c("Knowledge", "Behavior"), bold = TRUE) %>%
  align(j = c("T0_Med_IQR", "T1_Med_IQR"), align = "center", part = "all") %>%
  autofit() %>%
  flextable::save_as_docx(path = file.path(OUTPUT_DIR, "Changes_in_Knowledge_and_Adherence.docx"))
cat("Saved: Knowledge and behavior change table saved\n")


make_change_row <- function(df, var_T0, var_T1, label, is_total = FALSE, domain = NA) {
  valid_idx <- complete.cases(df[[var_T0]], df[[var_T1]])
  x0 <- df[[var_T0]][valid_idx]
  x1 <- df[[var_T1]][valid_idx]
  if (length(x0) == 0) return(NULL)

  wt <- wilcox.test(x0, x1, paired = TRUE)
  diffs <- x1 - x0
  n_total <- length(diffs)

  imp <- sprintf("%.1f%%", sum(diffs > 0) / n_total * 100)
  unc <- sprintf("%.1f%%", sum(diffs == 0) / n_total * 100)
  dec <- sprintf("%.1f%%", sum(diffs < 0) / n_total * 100)

  if (is_total) {
    t0_out <- sprintf("%.1f [%.1f, %.1f]", median(x0), quantile(x0, 0.25), quantile(x0, 0.75))
    t1_out <- sprintf("%.1f [%.1f, %.1f]", median(x1), quantile(x1, 0.25), quantile(x1, 0.75))
  } else {
    t0_out <- "-"
    t1_out <- "-"
  }

  data.frame(
    Variable       = label,
    Domain         = domain,
    Is_Total       = is_total,
    T0_Med_IQR     = t0_out,
    T1_Med_IQR     = t1_out,
    Improved       = imp,
    Unchanged      = unc,
    Declined       = dec,
    Wilcoxon_P_raw = wt$p.value,
    stringsAsFactors = FALSE
  )
}

var_list <- list(
  list("K_Total_T0", "K_Total_T1", "Knowledge total score (0-59)", TRUE,  "Knowledge"),
  list("K_Smoke_T0", "K_Smoke_T1", "  Smoking", FALSE, "Knowledge"),
  list("K_Alcohol_T0", "K_Alcohol_T1", "  Alcohol", FALSE, "Knowledge"),
  list("K_Obesity_T0", "K_Obesity_T1", "  Obesity", FALSE, "Knowledge"),
  list("K_PA_T0", "K_PA_T1", "  Physical activity", FALSE, "Knowledge"),
  list("K_Diet_T0", "K_Diet_T1", "  Diet", FALSE, "Knowledge"),
  list("K_Occu_T0", "K_Occu_T1", "  Occupation", FALSE, "Knowledge"),
  list("K_Infection_T0", "K_Infection_T1", "  Infection", FALSE, "Knowledge"),
  list("K_Sun_T0", "K_Sun_T1", "  Sun exposure", FALSE, "Knowledge"),
  list("K_Air_T0", "K_Air_T1", "  Air pollution", FALSE, "Knowledge"),
  list("K_Breast_T0", "K_Breast_T1", "  Breastfeeding", FALSE, "Knowledge"),
  list("total_score1_T0", "total_score1_T1", "WCRF behavior score (0-8)", TRUE, "Behavior"),
  list("smoke_score_T0", "smoke_score_T1", "  Smoking", FALSE, "Behavior"),
  list("BMI_score_T0", "BMI_score_T1", "  BMI", FALSE, "Behavior"),
  list("pa_score_T0", "pa_score_T1", "  Physical activity", FALSE, "Behavior"),
  list("fruvege_score_T0", "fruvege_score_T1", "  Fruit & vegetables", FALSE, "Behavior"),
  list("wholflour_score_T0", "wholflour_score_T1", "  Whole grains", FALSE, "Behavior"),
  list("fastfood_score_T0", "fastfood_score_T1", "  Fast food", FALSE, "Behavior"),
  list("redmeat_score_T0", "redmeat_score_T1", "  Red & processed meat", FALSE, "Behavior"),
  list("soda_score_T0", "soda_score_T1", "  Sugar-sweetened beverage", FALSE, "Behavior"),
  list("alcohol_score_T0", "alcohol_score_T1", "  Alcohol", FALSE, "Behavior")
)

rows <- lapply(var_list, function(v) {
  if (all(c(v[[1]], v[[2]]) %in% names(df))) {
    make_change_row(df, v[[1]], v[[2]], v[[3]], v[[4]], v[[5]])
  }
})
tbl2_raw <- do.call(rbind, Filter(Negate(is.null), rows))

tbl2_raw$Wilcoxon_P_holm <- tbl2_raw$Wilcoxon_P_raw

for (dom in c("Knowledge", "Behavior")) {
  idx <- which(tbl2_raw$Domain == dom & !tbl2_raw$Is_Total)
  tbl2_raw$Wilcoxon_P_holm[idx] <- p.adjust(tbl2_raw$Wilcoxon_P_raw[idx], method = "holm")
}

cat("\n===== Holm-adjusted p-value diagnostics=====\n")
print(
  tbl2_raw[, c("Domain", "Variable", "Wilcoxon_P_raw", "Wilcoxon_P_holm")] %>%
    mutate(
      Wilcoxon_P_raw  = round(Wilcoxon_P_raw, 4),
      Wilcoxon_P_holm = round(Wilcoxon_P_holm, 4),
      Sig_before = Wilcoxon_P_raw < 0.05,
      Sig_after  = Wilcoxon_P_holm < 0.05,
      Flip       = Sig_before != Sig_after
    ),
  row.names = FALSE
)

USE_HOLM_CORRECTION <- TRUE

tbl2 <- tbl2_raw
tbl2$Wilcoxon_P <- style_pvalue(
  if (USE_HOLM_CORRECTION) tbl2$Wilcoxon_P_holm else tbl2$Wilcoxon_P_raw,
  digits = 3
)

empty_row <- function(domain_name) {
  data.frame(Variable = domain_name, T0_Med_IQR = NA, T1_Med_IQR = NA,
             Improved = NA, Unchanged = NA, Declined = NA, Wilcoxon_P = NA)
}

tbl2_display <- bind_rows(
  empty_row("Knowledge"),
  tbl2 %>% filter(Domain == "Knowledge") %>% select(Variable, T0_Med_IQR, T1_Med_IQR, Improved, Unchanged, Declined, Wilcoxon_P),
  empty_row("Behavior"),
  tbl2 %>% filter(Domain == "Behavior") %>% select(Variable, T0_Med_IQR, T1_Med_IQR, Improved, Unchanged, Declined, Wilcoxon_P)
)

tbl2_display %>%
  flextable() %>%
  set_caption(caption = paste0(
    "Changes in cancer prevention knowledge and behaviors from baseline to follow-up.",
    if (USE_HOLM_CORRECTION) " P-values were Holm-adjusted within the knowledge (10 comparisons) and behavior (9 comparisons) subdomain families; total scores are unadjusted single tests." else ""
  )) %>%
  set_header_labels(
    Variable   = "Variable",
    T0_Med_IQR = "Baseline Median [IQR]",
    T1_Med_IQR = "Follow-up Median [IQR]",
    Improved   = "Improved",
    Unchanged  = "Unchanged",
    Declined   = "Declined",
    Wilcoxon_P = if (USE_HOLM_CORRECTION) "Wilcoxon p (Holm)" else "Wilcoxon p"
  ) %>%
  bold(i = ~ Variable %in% c("Knowledge", "Behavior"), bold = TRUE) %>%
  align(j = c("T0_Med_IQR", "T1_Med_IQR"), align = "center", part = "all") %>%
  autofit() %>%
  flextable::save_as_docx(path = file.path(OUTPUT_DIR, "Changes_in_Knowledge_and_Adherence.docx"))

cat("Saved: Change summaries saved (", if (USE_HOLM_CORRECTION) "Holm-adjusted" else "unadjusted", ")\n")

saveRDS(tbl2, file.path(OUTPUT_DIR, "Knowledge_Behavior_Changes.rds"))


library(moments)
library(dplyr)
library(flextable)
library(officer)

ceiling_stats <- function(df, var_T0, label, max_score) {
  x <- df[[var_T0]]
  x <- x[!is.na(x)]
  n <- length(x)
  data.frame(
    Domain        = label,
    Max_Score     = max_score,
    Median_IQR    = sprintf("%.1f [%.1f, %.1f]", median(x), quantile(x, 0.25), quantile(x, 0.75)),
    Prop_MaxScore = sprintf("%.1f%%", sum(x == max_score) / n * 100),
    Skewness      = round(skewness(x), 3),
    stringsAsFactors = FALSE
  )
}


regression_to_mean <- function(df, var_T0, var_T1, label, max_score, ceiling_threshold = 0.55) {
  valid_idx <- complete.cases(df[[var_T0]], df[[var_T1]])
  x0    <- df[[var_T0]][valid_idx]
  x1    <- df[[var_T1]][valid_idx]
  diffs <- x1 - x0
  n     <- length(x0)

  prop_max <- sum(x0 == max_score) / n
  if (prop_max >= ceiling_threshold) {
    message(sprintf("⚠ Skipping '%s': %.1f%% at max score (ceiling too severe for stratification)",
                    label, prop_max * 100))
    return(NULL)
  }

  cuts <- unique(quantile(x0, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE))

  if (length(cuts) < 4) {
    cuts  <- unique(quantile(x0, probs = c(0, 0.5, 1), na.rm = TRUE))
    labs  <- c("Low", "High")[seq_len(length(cuts) - 1)]
  } else {
    labs  <- c("Low", "Medium", "High")
  }

  groups <- cut(x0, breaks = cuts, include.lowest = TRUE, labels = labs)

  data.frame(x0, x1, diffs, groups) %>%
    group_by(groups) %>%
    summarise(
      N         = n(),
      T0_Median = round(median(x0), 2),
      T1_Median = round(median(x1), 2),
      Improved  = sprintf("%.1f%%", sum(diffs > 0) / n() * 100),
      Unchanged = sprintf("%.1f%%", sum(diffs == 0) / n() * 100),
      Declined  = sprintf("%.1f%%", sum(diffs < 0) / n() * 100),
      .groups   = "drop"
    ) %>%
    mutate(Domain = label) %>%
    rename(Stratum = groups) %>%
    select(Domain, Stratum, N, T0_Median, T1_Median, Improved, Unchanged, Declined)
}


domain_max_scores <- list(
  list("K_Smoke_T0",     "Smoking",           4),
  list("K_Alcohol_T0",   "Alcohol",           2),
  list("K_Obesity_T0",   "Obesity",           5),
  list("K_PA_T0",        "Physical activity", 2),
  list("K_Diet_T0",      "Diet",             16),
  list("K_Occu_T0",      "Occupation",        5),
  list("K_Infection_T0", "Infection",         9),
  list("K_Sun_T0",       "Sun exposure",      6),
  list("K_Air_T0",       "Air pollution",     7),
  list("K_Breast_T0",    "Breastfeeding",     3),
  list("K_Total_T0",     "Knowledge Total",  59)
)


rtm_domains <- list(
  list("K_Smoke_T0",     "K_Smoke_T1",     "Smoking",           4),
  list("K_Alcohol_T0",   "K_Alcohol_T1",   "Alcohol",           2),
  list("K_Obesity_T0",   "K_Obesity_T1",   "Obesity",           5),
  list("K_PA_T0",        "K_PA_T1",        "Physical activity", 2),
  list("K_Diet_T0",      "K_Diet_T1",      "Diet",             16),
  list("K_Occu_T0",      "K_Occu_T1",      "Occupation",        5),
  list("K_Infection_T0", "K_Infection_T1", "Infection",         9),
  list("K_Sun_T0",       "K_Sun_T1",       "Sun exposure",      6),
  list("K_Air_T0",       "K_Air_T1",       "Air pollution",     7),
  list("K_Breast_T0",    "K_Breast_T1",    "Breastfeeding",     3)
)

rtm_rows <- lapply(rtm_domains, function(v) {
  if (all(c(v[[1]], v[[2]]) %in% names(df))) {
    regression_to_mean(df, v[[1]], v[[2]], v[[3]], v[[4]])
  }
})
tbl_rtm <- do.call(rbind, Filter(Negate(is.null), rtm_rows))


tbl_ceiling <- do.call(rbind, Filter(Negate(is.null), lapply(domain_max_scores, function(v) {
  if (v[[1]] %in% names(df)) ceiling_stats(df, v[[1]], v[[2]], v[[3]])
})))

tbl_rtm <- do.call(rbind, Filter(Negate(is.null), lapply(rtm_domains, function(v) {
  if (all(c(v[[1]], v[[2]]) %in% names(df))) {
    regression_to_mean(df, v[[1]], v[[2]], v[[3]], v[[4]])
  }
})))

ft_ceiling <- tbl_ceiling %>%
  flextable() %>%
  set_header_labels(
    Domain        = "Domain",
    Max_Score     = "Max Score",
    Median_IQR    = "Median [IQR]",
    Prop_MaxScore = "% at Max Score",
    Skewness      = "Skewness (T0)"
  ) %>%
  bold(i = ~ Domain == "Knowledge Total") %>%
  bg(i = ~ as.numeric(gsub("%", "", Prop_MaxScore)) >= 15,
     j = "Prop_MaxScore", bg = "#FFE0B2") %>%
  align(j = c("Max_Score","Median_IQR","Prop_MaxScore","Skewness"),
        align = "center", part = "all") %>%
  bold(part = "header") %>%
  bg(part = "header", bg = "#D9E1F2") %>%
  autofit()

ft_rtm <- tbl_rtm %>%
  flextable() %>%
  set_header_labels(
    Domain    = "Domain",
    Stratum   = "Baseline Stratum",
    N         = "N",
    T0_Median = "T0 Median",
    T1_Median = "T1 Median",
    Improved  = "Improved",
    Unchanged = "Unchanged",
    Declined  = "Declined"
  ) %>%
  merge_v(j = "Domain") %>%
  valign(j = "Domain", valign = "top") %>%
  align(j = c("N","T0_Median","T1_Median","Improved","Unchanged","Declined"),
        align = "center", part = "all") %>%
  bold(part = "header") %>%
  bg(part = "header", bg = "#D9E1F2") %>%
  autofit()

out_path <- file.path(OUTPUT_DIR, "Knowledge_Ceiling_and_Baseline_Strata.docx")

read_docx() %>%
  body_add_par("Ceiling effect assessment of knowledge domains at baseline (T0).",
               style = "heading 2") %>%
  body_add_flextable(ft_ceiling) %>%
  body_add_par(paste(
    "Note: % at Max Score = proportion of participants achieving the maximum possible score.",
    "Negative skewness indicates scores clustered near the maximum (ceiling).",
    "A ceiling effect is conventionally flagged when % at Max Score ≥ 15–20%."
  ), style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("Changes in knowledge scores by baseline score level.",
               style = "heading 2") %>%
  body_add_flextable(ft_rtm) %>%
  body_add_par(paste(
    "Note: Participants were divided into tertiles (Low/Medium/High) by T0 baseline score.",
    "Domains with a ceiling effect too severe for meaningful stratification",
    "(% at Max Score ≥ 55%: Smoking, Alcohol, Obesity, Physical activity, Occupation, Sun exposure)",
    "are excluded from this table; their distributional characteristics are reported in the ceiling effect assessment",
    "Improved = T1 > T0; Unchanged = T1 = T0; Declined = T1 < T0.",
    "Changes within baseline strata cannot distinguish regression to the mean from other sources of change."
  ), style = "Normal") %>%
  print(target = out_path)

cat("Saved: Knowledge ceiling and baseline-strata summaries saved to:", out_path, "\n")


base_all_scored <- base_raw %>%
  score_knowledge() %>%
  filter(!is.na(K_Total), !is.na(total_score1))

excluded_group <- base_all_scored %>%
  filter(!user_id %in% df$user_id) %>%
  mutate(Group = "Excluded")

included_group <- df %>%
  select(user_id, age, sex, bs_edu, familyIncome, SES_3 = SES_3_T0,
         oe_m_current = oe_m_current_T0, mh_chronic = mh_chronic_T0,
         total_score1 = total_score1_T0, K_Total = K_Total_T0) %>%
  mutate(Group = "Included")

attrition_df <- bind_rows(included_group, excluded_group) %>%
  mutate(
    Group      = factor(Group, levels = c("Included", "Excluded")),
    sex_f      = factor(sex, levels = c(1, 0), labels = c("Male", "Female")),
    chronic_f  = factor(mh_chronic, levels = c(0, 1), labels = c("No", "Yes")),
    employed_f = factor(oe_m_current, levels = c(1, 0), labels = c("Employed", "Other")),
    bs_edu_f   = factor(bs_edu, levels = 0:3, labels = c("Primary", "Secondary", "Bachelor", "Postgrad")),
    inc_f      = factor(familyIncome, levels = 1:6, labels = c("<30k", "30-50k", "50-100k", "100-200k", "200-500k", ">500k"))
  )

format_pval_2sig <- function(x) {
  sapply(x, function(p) {
    if (is.na(p)) return(NA_character_)
    if (p < 0.001) return("<0.001")
    formatC(p, digits = 3, format = "f")
  })
}

set.seed(2026)
tbl_s1 <- attrition_df %>%
  select(Group, age, sex_f, bs_edu_f, inc_f, SES_3, chronic_f, employed_f, K_Total, total_score1) %>%
  tbl_summary(
    by = Group,
    type = list(
      c(age, SES_3, K_Total, total_score1) ~ "continuous",
      c(sex_f, chronic_f, employed_f) ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{median} [{p25}, {p75}]",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_categorical() ~ c(0, 1),
      all_continuous() ~ 1
    ),
    label = list(
      age ~ "Age (years)",
      sex_f ~ "Sex",
      bs_edu_f ~ "Education",
      inc_f ~ "Annual income (CNY)",
      SES_3 ~ "Subjective SES (1-10)",
      chronic_f ~ "Chronic disease",
      employed_f ~ "Currently employed",
      K_Total ~ "Knowledge total score (0-59)",
      total_score1 ~ "WCRF behavior score (0-8)"
    ),
    missing = "no"
  ) %>%
  add_overall() %>%
  add_p(
    test = list(all_continuous() ~ "wilcox.test", all_categorical() ~ "chisq.test"),
    test.args = all_categorical() ~ list(simulate.p.value = TRUE),
    pvalue_fun = format_pval_2sig
  ) %>%
  modify_caption("**Comparison of included and excluded participants**") %>%
  bold_labels()

tbl_s1 %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = file.path(OUTPUT_DIR, "Included_vs_Excluded_Comparison.docx"))

cat("Saved: Participant comparison table saved\n")

saveRDS(as_tibble(tbl1), file.path(OUTPUT_DIR, "Baseline_Characteristics.rds"))

saveRDS(tbl2, file.path(OUTPUT_DIR, "Knowledge_Behavior_Changes.rds"))

saveRDS(as_tibble(tbl_s1), file.path(OUTPUT_DIR, "Included_vs_Excluded_Comparison.rds"))

n_table_s1_total <- nrow(included_group) + nrow(excluded_group)
cat("Included/excluded cohort total N:", n_table_s1_total, "\n")
cat("Matches eligible baseline N:", n_table_s1_total == n_base_valid, "\n")


library(survey)
library(cobalt)
library(boot)
library(dplyr)
library(flextable)
library(officer)

attrition_df$included_bin <- as.numeric(attrition_df$Group == "Included")

ps_model <- glm(
  included_bin ~ age + sex_f + bs_edu_f + inc_f + SES_3 +
    chronic_f + employed_f + K_Total + total_score1,
  data = attrition_df,
  family = binomial
)

attrition_df$ps  <- predict(ps_model, type = "response")
attrition_df$ipw <- ifelse(attrition_df$Group == "Included",
                           1 / attrition_df$ps,
                           1)

cat("\n=== IPW weight distribution among included participants===\n")
print(summary(attrition_df$ipw[attrition_df$Group == "Included"]))

bal_out <- bal.tab(ps_model,
                   data = attrition_df,
                   weights = attrition_df$ipw,
                   method = "weighting",
                   un = TRUE)

bal_df <- as.data.frame(bal_out$Balance) %>%
  tibble::rownames_to_column("Covariate") %>%
  select(Covariate, Diff.Un, Diff.Adj) %>%
  rename(`SMD (Unweighted)` = Diff.Un, `SMD (Weighted)` = Diff.Adj) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

w_inc         <- attrition_df$ipw[attrition_df$Group == "Included"]
ess_unweighted <- length(w_inc)
ess_weighted   <- sum(w_inc)^2 / sum(w_inc^2)

cat("\nEffective sample size: unweighted =", ess_unweighted,
    ", weighted =", round(ess_weighted, 1), "\n")

ipw_lookup <- attrition_df %>%
  filter(Group == "Included") %>%
  select(user_id, ipw)

df_weighted <- cohort_clean_final %>%
  left_join(ipw_lookup, by = "user_id")

cat("Unmatched participants:", sum(is.na(df_weighted$ipw)), "\n")

diff_var  <- df_weighted$K_Total_T1 - df_weighted$K_Total_T0
w         <- df_weighted$ipw
boot_data <- data.frame(diff = diff_var, w = w)

set.seed(2026)

wtd_mean_boot <- function(data, indices) {
  d <- data[indices, ]
  weighted.mean(d$diff, w = d$w)
}
boot_out_weighted <- boot(boot_data, wtd_mean_boot, R = 2000)
ci_weighted <- boot.ci(boot_out_weighted, type = "perc")$percent[4:5]

mean_boot <- function(data, indices) {
  d <- data[indices, ]
  mean(d$diff)
}
boot_out_unweighted <- boot(boot_data, mean_boot, R = 2000)
ci_unweighted <- boot.ci(boot_out_unweighted, type = "perc")$percent[4:5]

panelB_df <- data.frame(
  Analysis = c("Unweighted", "IPW-weighted"),
  `Mean change (T1-T0)` = c(
    round(mean(diff_var, na.rm = TRUE), 2),
    round(weighted.mean(diff_var, w = w, na.rm = TRUE), 2)
  ),
  `95% CI (bootstrap, R=2000)` = c(
    sprintf("(%.2f, %.2f)", ci_unweighted[1], ci_unweighted[2]),
    sprintf("(%.2f, %.2f)", ci_weighted[1], ci_weighted[2])
  ),
  check.names = FALSE
)

ft_balance <- bal_df %>%
  flextable() %>%
  align(j = c("SMD (Unweighted)", "SMD (Weighted)"), align = "center", part = "all") %>%
  bold(part = "header") %>%
  bg(part = "header", bg = "#D9E1F2") %>%
  autofit()

ft_meanchange <- panelB_df %>%
  flextable() %>%
  align(j = c("Mean change (T1-T0)", "95% CI (bootstrap, R=2000)"), align = "center", part = "all") %>%
  bold(part = "header") %>%
  bg(part = "header", bg = "#D9E1F2") %>%
  autofit()

out_path_s4 <- file.path(OUTPUT_DIR, "IPW_Balance_and_Knowledge_Change.docx")

read_docx() %>%
  body_add_par("Covariate balance before and after inverse probability weighting (IPW).",
               style = "heading 2") %>%
  body_add_flextable(ft_balance) %>%
  body_add_par(paste(
    "Note: SMD = standardized mean difference. Propensity scores for retention were estimated",
    "via logistic regression on all measured baseline covariates."
  ), style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("Weighted versus unweighted change in knowledge total score (T1−T0).",
               style = "heading 2") %>%
  body_add_flextable(ft_meanchange) %>%
  body_add_par(paste(
    sprintf("Note: Because inverse probability weighting adjusts expectation-based estimators, weighted analyses are based on mean change rather than the median-based Wilcoxon test used in the primary analysis of within-person change. Effective sample size after weighting: %.1f (unweighted N=%d, a %.1f%% reduction).",
            ess_weighted, ess_unweighted, (1 - ess_weighted/ess_unweighted) * 100)
  ), style = "Normal") %>%
  print(target = out_path_s4)

cat("Saved: Supplementary IPW descriptive results saved to:", out_path_s4, "\n")


ipw_export <- attrition_df %>%
  filter(Group == "Included") %>%
  select(user_id, ps, ipw)

df_805_check <- readRDS(file.path(OUTPUT_DIR, "cohort_clean_N805.rds"))
n_match <- sum(df_805_check$user_id %in% ipw_export$user_id)
cat("Participants matched to IPW weights:", n_match, "\n")
if (n_match != nrow(df_805_check)) {
  warning("Some participants did not match the IPW weight file; verify study ID alignment.")
}

saveRDS(ipw_export, file.path(OUTPUT_DIR, "ipw_weights_N805.rds"))
cat("Saved: IPW weights saved for downstream scripts:",
    file.path(OUTPUT_DIR, "ipw_weights_N805.rds"), "\n")
