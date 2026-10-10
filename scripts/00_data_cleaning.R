# Data preparation and analytic cohort selection.
# No individual-level data are exported publicly.


library(dplyr)
library(haven)
library(lubridate)

DATA_DIR <- Sys.getenv("CEBP_DATA_DIR", unset = file.path(getwd(), "data"))
OUTPUT_DIR <- Sys.getenv("CEBP_OUTPUT_DIR", unset = file.path(getwd(), "output"))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


score_knowledge <- function(df) {
  df %>% mutate(
    K_Smoke      = sq1 + sq2 + sq18,
    K_Obesity    = sq3 + sq19 + sq20,
    K_PA         = sq4 + sq5,
    K_Alcohol    = sq6 + sq7,
    K_PlantFiber = sq12 + sq21,
    K_FriedMeat  = sq29,
    K_Sugar      = sq11,
    K_HotFood    = sq10 + sq28,
    K_Hazard     = sq8  + sq9  + sq27,
    K_Infection  = sq13 + sq17 + sq22 + sq26 + sq30,
    K_Occu       = sq14 + sq23 + sq24,
    K_Sun        = sq15 + sq25 + sq31,
    K_Air        = sq16 + sq32 + sq33,
    K_Breast     = sq34,
    K_Diet       = sq8+sq9+sq10+sq11+sq12+sq21+sq27+sq28+sq29,
    K_Total      = sq1 +sq2 +sq3 +sq4 +sq5 + sq6 +sq7 +sq8 +sq9 +sq10 +
      sq11 +sq12 +sq13 +sq14 +sq15 + sq16 +sq17 +sq18 +sq19 +sq20 +
      sq21 +sq22 +sq23 +sq24 +sq25 + sq26 +sq27 +sq28 +sq29 +sq30 +
      sq31 +sq32 +sq33 +sq34
  )
}

score_behavior <- function(df) {
  df %>% mutate(
    hotfood_score = case_when(
      hotfood == 0 ~ 1,
      hotfood == 1 ~ 0.5,
      hotfood >= 2 ~ 0,
      TRUE         ~ NA_real_
    ),
    dessert_score = case_when(
      dh_dessert == 0 ~ 1,
      dh_dessert == 1 ~ 0.5,
      dh_dessert >= 2 ~ 0,
      TRUE            ~ NA_real_
    ),
    fried_score = case_when(
      dh_dfry == 0 ~ 1,
      dh_dfry == 1 ~ 0.5,
      dh_dfry >= 2 ~ 0,
      TRUE         ~ NA_real_
    ),
    betel_score = ifelse(pm_betel == 0, 1, 0),
    hepb_score  = ifelse(hep_B_va == 1, 1, 0),
    hpv_score   = ifelse(hpv_va_pm == 1, 1, 0),

    wcrf7 = BMI_score + pa_score + fruvege_score + wholflour_score +
      fastfood_score + redmeat_score + soda_score + alcohol_score,
    wcrf8 = wcrf7 + smoke_score,
    wcrf9 = wcrf8 + hotfood_score
  )
}

cat("\n=== Data preparation ===\n")
required_inputs <- file.path(DATA_DIR, c("base.sas7bdat", "follow.sas7bdat"))
if (!all(file.exists(required_inputs))) {
  stop("Missing private SAS inputs: ", paste(basename(required_inputs[!file.exists(required_inputs)]), collapse = ", "),
       ". Place authorized files in data/ or set CEBP_DATA_DIR; never commit participant data.")
}
base_raw   <- read_sas(file.path(DATA_DIR, "base.sas7bdat"))
follow_raw <- read_sas(file.path(DATA_DIR, "follow.sas7bdat"))

f1_target_names <- follow_raw %>% select(ends_with("_f1")) %>% names() %>% gsub("_f1$", "", .)
conflict_vars <- intersect(names(follow_raw), f1_target_names)
follow_clean <- follow_raw %>%
  select(-all_of(conflict_vars)) %>%
  select(user_id, kn_date_start_f1, kn_date_end_f1, total_score1_f1,
         starts_with("sq"), starts_with("q")) %>%
  rename_with(~ gsub("_f1", "", .), everything())

raw_q_base <- base_raw %>%
  select(user_id, all_of(paste0("q", 1:26)))

raw_q_follow <- follow_raw %>%
  select(user_id, all_of(paste0("q", 1:26, "_f1"))) %>%
  rename_with(~ gsub("_f1", "", .), -user_id)

base_scored   <- base_raw %>% score_knowledge() %>% score_behavior()
follow_scored <- follow_clean %>% score_knowledge()

covariates <- c("user_id", "age", "sex", "bs_edu", "familyIncome", "bmi", "affiliation_grp")

dat_T0 <- base_scored %>% rename_with(~ paste0(., "_T0"), .cols = -all_of(covariates))
dat_T1 <- follow_scored %>% rename_with(~ paste0(., "_T1"), .cols = -user_id)

merged_split <- inner_join(dat_T0, dat_T1, by = "user_id")


q_cols <- paste0("q", 1:26)

raw_q_base_num <- raw_q_base %>%
  mutate(across(all_of(q_cols), ~ suppressWarnings(as.numeric(.))))

raw_q_follow_num <- raw_q_follow %>%
  mutate(across(all_of(q_cols), ~ suppressWarnings(as.numeric(.))))

sd_T0 <- raw_q_base_num %>%
  mutate(q_sd_T0 = apply(select(., all_of(q_cols)), 1, sd, na.rm = TRUE)) %>%
  select(user_id, q_sd_T0)

sd_T1 <- raw_q_follow_num %>%
  mutate(q_sd_T1 = apply(select(., all_of(q_cols)), 1, sd, na.rm = TRUE)) %>%
  select(user_id, q_sd_T1)

cat("\n=== Straight-lining diagnostics based on responses q1–q26===\n")
cat("T0 zero-variance responses:", sum(sd_T0$q_sd_T0 == 0, na.rm=TRUE), "participants\n")
cat("T1 zero-variance responses:", sum(sd_T1$q_sd_T1 == 0, na.rm=TRUE), "participants\n")
cat("T0 low-variance responses (SD <= 0.1):", sum(sd_T0$q_sd_T0 <= 0.1, na.rm=TRUE), "participants\n")
cat("T1 low-variance responses (SD <= 0.1):", sum(sd_T1$q_sd_T1 <= 0.1, na.rm=TRUE), "participants\n")

df_filtered <- merged_split %>%
  left_join(sd_T0, by = "user_id") %>%
  left_join(sd_T1, by = "user_id") %>%
  mutate(
    t0_clean  = na_if(kn_date_start_T0, ""),
    t1_clean  = na_if(kn_date_start_T1, ""),
    Date_T0   = as.POSIXct(t0_clean, format = "%Y-%m-%d %H:%M:%S"),
    Date_T1   = as.POSIXct(t1_clean, format = "%Y-%m-%d %H:%M:%S"),
    Days_Diff = as.numeric(difftime(Date_T1, Date_T0, units = "days"))
  ) %>%
  filter(
    Days_Diff >= 30,
    !is.na(K_Total_T0), !is.na(K_Total_T1),
    !is.na(total_score1_T0), !is.na(total_score1_T1),
    is.na(q_sd_T0) | q_sd_T0 > 0,
    is.na(q_sd_T1) | q_sd_T1 > 0
  )

cat("\nN after minimum follow-up and response-quality exclusions =", nrow(df_filtered), "\n")

follow_behavior_T1 <- follow_raw %>%
  select(
    user_id,
    hotfood = hotfood_f1, dh_dessert = dh_dessert_f1, dh_dfry = dh_dfry_f1,
    pm_betel = pm_betel_f1, hep_B_va = hep_B_va_f1, hpv_va_pm = hpv_va_pm_f1,
    BMI_score = BMI_score_f1, pa_score = pa_score_f1, fruvege_score = fruvege_score_f1,
    wholflour_score = wholflour_score_f1, fastfood_score = fastfood_score_f1,
    redmeat_score = redmeat_score_f1, soda_score = soda_score_f1,
    smoke_score = smoke_score_f1, alcohol_score = alcohol_score_f1
  ) %>%
  score_behavior() %>%
  select(
    user_id, BMI_score, pa_score, fastfood_score, redmeat_score, soda_score,
    smoke_score, alcohol_score, hotfood_score, dessert_score, fried_score,
    betel_score, hepb_score, hpv_score, fruvege_score, wholflour_score,
    wcrf7, wcrf8, wcrf9
  ) %>%
  rename_with(~ paste0(., "_T1"), .cols = -user_id)

df_split <- df_filtered %>% left_join(follow_behavior_T1, by = "user_id")

cat("\n=== Missing data diagnostics ===\n")
required_vars <- c("K_Total_T0", "total_score1_T0", "K_Total_T1", "wcrf8_T1",
                   "age", "sex", "bs_edu", "familyIncome", "SES_3_T0", "mh_chronic_T0", "oe_m_current_T0")
missing_summary <- colSums(is.na(df_split %>% select(all_of(required_vars))))
print(missing_summary[missing_summary > 0])

cohort_clean_final <- df_split %>% filter(complete.cases(select(., all_of(required_vars))))
cat("\nFinal analytic sample: N =", nrow(cohort_clean_final))

saveRDS(cohort_clean_final, file.path(OUTPUT_DIR, "cohort_clean_N805.rds"))
cat(sprintf("Cleaned cohort saved to %s\n", OUTPUT_DIR))


cor_wcrf <- cor(cohort_clean_final$wcrf8_T0, cohort_clean_final$total_score1_T0, use = "complete.obs")
cat("\n[1] Compare R adherence score with SAS source score:\n")
cat("  Correlation r =", round(cor_wcrf, 5), "\n")
if(cor_wcrf < 0.999) {
  cat("  WARNING: source scores differ from recomputed scores:\n")
  message("Number of participants with scoring discrepancies: ",
          sum(abs(cohort_clean_final$wcrf8_T0 - cohort_clean_final$total_score1_T0) > 0.01, na.rm = TRUE))
  cat("  -> Check upstream SAS scoring and component definitions.\n")
} else {
  cat("  R and SAS adherence scores agree (r >= 0.999).\n")
}

cat("\nFinal analytic sample: N =", nrow(cohort_clean_final), "\n")


n_base_raw <- nrow(base_raw)

base_all_scored <- base_raw %>%
  score_knowledge() %>%
  filter(!is.na(K_Total), !is.na(total_score1))
n_base_valid   <- nrow(base_all_scored)
n_excl_base_incomplete <- n_base_raw - n_base_valid

n_follow_raw <- nrow(follow_raw)

n_matched_valid <- length(intersect(base_all_scored$user_id, follow_raw$user_id))
n_lost_to_followup <- n_base_valid - n_matched_valid

cat("Eligible baseline participants:", n_base_valid, "\n")
cat("Participants with follow-up:", n_matched_valid, "\n")
cat("Participants lost to follow-up:", n_lost_to_followup, "\n")
