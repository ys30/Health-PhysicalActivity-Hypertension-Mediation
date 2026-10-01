# NHANES 2017–March 2020 pre-pandemic
# Physical activity, waist circumference, and hypertension
#
# X: MVPA_me  ->  M: waist  ->  Y: HTN
# Covariates: age, sex, race/ethnicity, PIR, sedentary time, smoking, alcohol
# Adults: 18–62
#
# Note: this is an unweighted exploratory analysis of the selected complete-case
# sample. The model sequence is pathway/mediation-style, not a formal causal
# mediation estimator.

library(tidyverse)
library(haven)
library(DiagrammeR)
library(broom)
library(lmtest)
library(car)
library(pROC)
library(ResourceSelection)

# 1. Download public NHANES components ---------------------------------------

urls <- c(
  demo = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_DEMO.XPT",
  bmx  = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_BMX.XPT",
  bpx  = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_BPXO.XPT",
  paq  = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_PAQ.XPT",
  bpq  = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_BPQ.XPT",
  smq  = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_SMQ.XPT",
  alq  = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_ALQ.XPT"
)

read_nhanes <- function(url) {
  tf <- tempfile(fileext = ".xpt")
  download.file(url, tf, mode = "wb", quiet = TRUE)
  read_xpt(tf)
}

nhanes <- lapply(urls, read_nhanes)
demo <- nhanes$demo; bmx <- nhanes$bmx; bpx <- nhanes$bpx
paq <- nhanes$paq; bpq <- nhanes$bpq; smq <- nhanes$smq; alq <- nhanes$alq

# 2. Physical activity --------------------------------------------------------

paq_mvpa <- paq %>%
  transmute(
    SEQN,
    work_vig_days = if_else(PAQ605 == 1, as.numeric(PAQ610), 0),
    work_vig_min  = if_else(PAQ605 == 1, as.numeric(PAD615), 0),
    work_mod_days = if_else(PAQ620 == 1, as.numeric(PAQ625), 0),
    work_mod_min  = if_else(PAQ620 == 1, as.numeric(PAD630), 0),
    trans_days    = if_else(PAQ635 == 1, as.numeric(PAQ640), 0),
    trans_min     = if_else(PAQ635 == 1, as.numeric(PAD645), 0),
    rec_vig_days  = if_else(PAQ650 == 1, as.numeric(PAQ655), 0),
    rec_vig_min   = if_else(PAQ650 == 1, as.numeric(PAD660), 0),
    rec_mod_days  = if_else(PAQ665 == 1, as.numeric(PAQ670), 0),
    rec_mod_min   = if_else(PAQ665 == 1, as.numeric(PAD675), 0),
    sedentary_min_day = as.numeric(PAD680)
  ) %>%
  mutate(
    vig_min_week = work_vig_days * work_vig_min + rec_vig_days * rec_vig_min,
    mod_min_week = work_mod_days * work_mod_min +
      rec_mod_days * rec_mod_min + trans_days * trans_min,
    MVPA_me = mod_min_week + 2 * vig_min_week
  ) %>%
  select(SEQN, MVPA_me, sedentary_min_day)

# 3. Blood pressure outcome ---------------------------------------------------

bp <- bpx %>%
  transmute(
    SEQN,
    SBP = rowMeans(select(., BPXOSY1, BPXOSY2, BPXOSY3), na.rm = TRUE),
    DBP = rowMeans(select(., BPXODI1, BPXODI2, BPXODI3), na.rm = TRUE)
  ) %>%
  mutate(
    HTN = if_else(
      !is.na(SBP) & !is.na(DBP),
      if_else(SBP >= 130 | DBP >= 80, 1, 0),
      NA_real_
    )
  )

# 4. Mediator and covariates --------------------------------------------------

waist <- bmx %>% transmute(SEQN, waist = as.numeric(BMXWAIST))

smoking <- smq %>%
  transmute(SEQN, smoker = case_when(
    SMQ020 == 1 ~ 1,
    SMQ020 == 2 ~ 0,
    TRUE ~ NA_real_
  ))

alcohol_use <- alq %>%
  transmute(SEQN, alcohol = case_when(
    ALQ111 == 1 ~ 1,
    ALQ111 == 2 ~ 0,
    TRUE ~ NA_real_
  ))

bp_meds <- bpq %>%
  transmute(SEQN, bp_meds = case_when(
    BPQ050A == 1 ~ 1,
    BPQ050A == 2 ~ 0,
    TRUE ~ NA_real_
  ))

demographics <- demo %>%
  filter(between(RIDAGEYR, 18, 62)) %>%
  transmute(
    SEQN,
    age = as.numeric(RIDAGEYR),
    sex = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female")),
    race = factor(RIDRETH3),
    pir = as.numeric(INDFMPIR)
  )

# 5. Analysis dataset ---------------------------------------------------------

df <- demographics %>%
  left_join(paq_mvpa, by = "SEQN") %>%
  left_join(waist, by = "SEQN") %>%
  left_join(bp, by = "SEQN") %>%
  left_join(smoking, by = "SEQN") %>%
  left_join(alcohol_use, by = "SEQN") %>%
  left_join(bp_meds, by = "SEQN") %>%
  filter(
    !is.na(MVPA_me), is.finite(MVPA_me),
    !is.na(waist), is.finite(waist),
    !is.na(SBP), is.finite(SBP),
    !is.na(DBP), is.finite(DBP),
    !is.na(HTN),
    !is.na(sedentary_min_day), is.finite(sedentary_min_day),
    !is.na(smoker), !is.na(alcohol),
    !is.na(age), !is.na(sex), !is.na(race), !is.na(pir),
    !is.na(bp_meds)
  ) %>%
  filter(bp_meds == 0)

cat("Analytic n (no BP medications) =", nrow(df), "\n")

# 6. Descriptive outputs ------------------------------------------------------

sample_summary <- df %>%
  summarise(
    n = n(),
    age_mean = mean(age), age_sd = sd(age),
    pir_mean = mean(pir), pir_sd = sd(pir),
    MVPA_mean = mean(MVPA_me), MVPA_sd = sd(MVPA_me),
    waist_mean = mean(waist), waist_sd = sd(waist),
    SBP_mean = mean(SBP), SBP_sd = sd(SBP),
    HTN_prevalence = mean(HTN),
    sedentary_min_day_mean = mean(sedentary_min_day)
  )

dir.create("results", showWarnings = FALSE)
write.csv(sample_summary, "results/sample_summary_generated.csv", row.names = FALSE)

# 7. Pathway-oriented regression models --------------------------------------

model_M <- lm(
  waist ~ MVPA_me + sedentary_min_day + smoker + alcohol +
    age + sex + race + pir,
  data = df
)

model_H1 <- glm(
  HTN ~ MVPA_me + sedentary_min_day + smoker + alcohol +
    age + sex + race + pir,
  data = df,
  family = binomial()
)

model_H2 <- glm(
  HTN ~ MVPA_me + waist + sedentary_min_day + smoker + alcohol +
    age + sex + race + pir,
  data = df,
  family = binomial()
)

print(summary(model_M))
print(summary(model_H1))
print(summary(model_H2))

# Compare H1 with H2 and inspect change in MVPA coefficient.
get_mvpa_beta <- function(mod) {
  tidy(mod) %>% filter(term == "MVPA_me") %>% pull(estimate)
}

attenuation_table <- tibble(
  model = c("H1: without waist", "H2: with waist"),
  MVPA_beta = c(get_mvpa_beta(model_H1), get_mvpa_beta(model_H2)),
  AIC = c(model_H1$aic, model_H2$aic)
)

write.csv(attenuation_table, "results/mvpa_attenuation_generated.csv", row.names = FALSE)
print(attenuation_table)
print(anova(model_H1, model_H2, test = "Chisq"))

# 8. Diagnostics --------------------------------------------------------------

cat("\nVIF: waist model\n")
print(car::vif(model_M))
cat("\nBreusch-Pagan: waist model\n")
print(lmtest::bptest(model_M))

logit_diagnostics <- function(mod, name) {
  p_hat <- predict(mod, type = "response")
  y_obs <- mod$y
  roc_obj <- pROC::roc(y_obs, p_hat, quiet = TRUE)

  cat("\n", name, " AUC: ", as.numeric(pROC::auc(roc_obj)), "\n", sep = "")
  print(ResourceSelection::hoslem.test(y_obs, p_hat, g = 10))

  cook <- cooks.distance(mod)
  cutoff <- 4 / nrow(model.frame(mod))
  cat("Influential observations (Cook's D > 4/n): ",
      sum(cook > cutoff), "\n", sep = "")

  invisible(list(auc = as.numeric(pROC::auc(roc_obj)), cooks_d = cook))
}

diag_H1 <- logit_diagnostics(model_H1, "H1")
diag_H2 <- logit_diagnostics(model_H2, "H2")

# 9. Conceptual pathway diagram ----------------------------------------------

DiagrammeR::grViz("
digraph dag {
  graph [layout = dot, rankdir = LR]
  MVPA [label='MVPA (X)']
  W [label='Waist circumference (M)']
  HTN [label='Hypertension (Y)']
  C [label='Age, sex, race/ethnicity, PIR']
  L [label='Sedentary time, smoking, alcohol']
  MVPA -> W
  W -> HTN
  MVPA -> HTN
  C -> MVPA
  C -> W
  C -> HTN
  L -> MVPA
  L -> W
  L -> HTN
}
")
