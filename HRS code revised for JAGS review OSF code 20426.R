# ==============================================================================
# Stroke Mortality Analysis: Comparing 1996-2006 vs 2008-2018 periods

# ==============================================================================
options(scipen=999)

rm(list=ls())
library(haven)        # to read stata data
library(tidyverse)    # data manipulation
library(dtms)         # discrete time multistate models
library(msm)          # continuous multistate model
library(elect)        # continuous multistate model
library(survival)     # survival analysis
library(survminer)    # survival analysis
library(tableone)     # for Table 1

HRSdata <- read_dta("path/randhrs1992_2020v1.dta")

# ------------------------------------------------------------------------------
# 1. DATA PREPARATION -
# ------------------------------------------------------------------------------
 
# Select relevant variables from HRS data
# 
 deathdata <- HRSdata %>%
   select(
     # ID and demographics
     hhidpn, ragender, rabyear, rabmonth, raeduc, raedyrs,
     # Race/ethnicity (time-invariant)
     raracem, rahispan,
     
     # Death information
     radmonth, radyear, raddate,
     
     # Health conditions (waves 2-15)
     starts_with("r") & (ends_with("stroke") | ends_with("strok") |
                           ends_with("stroks") | ends_with("heart") |
                           ends_with("hibpe") | ends_with("diabe") |
                           ends_with("mobil") | ends_with("smokev") |
                           ends_with("mstat") | ends_with("iwstat") |
                           ends_with("bmi") | ends_with("higov") |
                           ends_with("prpcnt") | ends_with("covr"))
 )
  
# ------------------------------------------------------------------------------
# 2. RESHAPE TO LONG FORMAT
# ------------------------------------------------------------------------------

deathdatasetlong <- deathdata %>%
  pivot_longer(
    cols = matches("r\\d+(stroks|heart|strok|stroke|hibpe|iwstat|diabe|mobil|mstat|smokev|bmi|higov|prpcnt|covr)"),
    names_to = c("wave", ".value"),
    names_pattern = "r(\\d+)(stroks|heart|strok|stroke|hibpe|iwstat|diabe|mobil|mstat|smokev|bmi|higov|prpcnt|covr)"
  ) %>%
  mutate(
    wave = recode(wave,
                  `2` = 1994, `3` = 1996, `4` = 1998, `5` = 2000,
                  `6` = 2002, `7` = 2004, `8` = 2006, `9` = 2008,
                  `10` = 2010, `11` = 2012, `12` = 2014, `13` = 2016,
                  `14` = 2018, `15` = 2020),
    # Clean death dates
    radyear = if_else(radyear > 2020, NA_real_, radyear),
    raddate = if_else(is.na(radyear), NA_real_, raddate)
  ) %>%
  # Remove rows with no health information
  filter(!(is.na(heart) & is.na(strok) & is.na(stroks) &
             is.na(diabe) & is.na(smokev) & is.na(mstat)))

# ------------------------------------------------------------------------------
# 3. CALCULATE AGE AND SURVIVAL TIME
# ------------------------------------------------------------------------------

deathdatasetlong <- deathdatasetlong %>%
  group_by(hhidpn) %>%
  filter(radyear > 1993 | is.na(radyear)) %>%
  mutate(
    first_wave = min(wave, na.rm = TRUE),
    
    # Month fraction for precise age calculation
    month_fractionb = case_when(
      rabmonth == 1 ~ 15/365, rabmonth == 2 ~ 45/365, rabmonth == 3 ~ 74/365,
      rabmonth == 4 ~ 105/365, rabmonth == 5 ~ 135/365, rabmonth == 6 ~ 166/365,
      rabmonth == 7 ~ 196/365, rabmonth == 8 ~ 227/365, rabmonth == 9 ~ 258/365,
      rabmonth == 10 ~ 288/365, rabmonth == 11 ~ 319/365, rabmonth == 12 ~ 349/365,
      TRUE ~ 0
    ),
    
    month_fractiond = case_when(
      radmonth == 1 ~ 15/365, radmonth == 2 ~ 45/365, radmonth == 3 ~ 74/365,
      radmonth == 4 ~ 105/365, radmonth == 5 ~ 135/365, radmonth == 6 ~ 166/365,
      radmonth == 7 ~ 196/365, radmonth == 8 ~ 227/365, radmonth == 9 ~ 258/365,
      radmonth == 10 ~ 288/365, radmonth == 11 ~ 319/365, radmonth == 12 ~ 349/365,
      TRUE ~ 0
    ),
    
    # Age calculations
    Age = wave - rabyear - month_fractionb,
    Agedeath = radyear - rabyear + month_fractiond,
    Agefirstwave = first_wave - rabyear - month_fractionb,
    
    # Survival time
    Survivalyears = if_else(
      (Agedeath - Age) < 0,
      Agedeath - Agefirstwave,
      Agedeath - Age
    ),
    Survivalmonth = round(Survivalyears * 12),
    # Filter for age 50+
    Age = ifelse(Age > 49 & (raddate > 0 | is.na(raddate)), Age, NA)
  ) %>%
  ungroup() %>%
  filter(!is.na(Age))

# ------------------------------------------------------------------------------
# 4. CREATE EVENT INDICATOR AND RECODE VARIABLES
# ------------------------------------------------------------------------------

deathdatasetlong <- deathdatasetlong %>%
  group_by(hhidpn) %>%
  mutate(event = ifelse(wave == max(wave) & Survivalyears > 0, 1, 0)) %>%
  ungroup() %>%
  mutate(
    event = replace_na(event, 0),
    # Recode binary variables
    strok = ifelse(strok == 1, 1, 0),
    smokev = ifelse(smokev == 1, 1, 0),
    hibpe = ifelse(hibpe == 1, 1, 0),
    diabe = ifelse(diabe == 1, 1, 0),
    heart = ifelse(heart == 1, 1, 0),
    
    # Create race/ethnicity variable
    race_eth = case_when(
      rahispan == 1 ~ "Hispanic",
      raracem == 1 ~ "White",
      raracem == 2 ~ "Black",
      raracem == 3 ~ "Other",
      TRUE ~ NA_character_
    ),
    race_eth = factor(race_eth, levels = c("White", "Black", "Hispanic", "Other")),
    
    # Create BMI - harmonize across wave formats
    # RAND HRS stores BMI differently across waves:
    #   Earlier waves: BMI in correct scale (e.g., 25 = BMI 25.0)
    #   Later waves: BMI * 10 (e.g., 250 = BMI 25.0)
    bmi_real = case_when(
      bmi > 100 ~ bmi / 10,
      bmi <= 100 ~ bmi,
      TRUE ~ NA_real_
    ),
    bmi_cat = case_when(
      bmi_real < 18.5 ~ "Underweight",
      bmi_real >= 18.5 & bmi_real < 25 ~ "Normal",
      bmi_real >= 25 & bmi_real < 30 ~ "Overweight",
      bmi_real >= 30 ~ "Obese",
      TRUE ~ NA_character_
    ),
    bmi_cat = factor(bmi_cat, levels = c("Normal", "Underweight", "Overweight", "Obese")),
    
    # NEW: Create insurance variable
    insurance_type = case_when(
      higov == 1 & prpcnt > 0  ~ "Both",
      higov == 1 & (prpcnt == 0 | is.na(prpcnt)) ~ "Government only",
      (higov == 0 | is.na(higov)) & prpcnt > 0  ~ "Private only",
      (higov == 0 | is.na(higov)) & (prpcnt == 0 | is.na(prpcnt)) ~ "Uninsured",
      TRUE ~ NA_character_
    ),
    insurance_type = factor(insurance_type,
                            levels = c("Private only", "Government only", "Both", "Uninsured"))
  ) %>%
  distinct(hhidpn, wave, .keep_all = TRUE) %>%
  mutate(
    strok = factor(strok, levels = c(0, 1), labels = c("No Stroke", "Stroke")),
    ragender = factor(ragender, levels = c(1, 2), labels = c("Male", "Female"))
  )

# Check what bmi actually looks like
summary(deathdatasetlong$bmi)

# See the distribution
quantile(deathdatasetlong$bmi, probs = seq(0, 1, 0.1), na.rm = TRUE)

# Check a handful of actual values
head(deathdatasetlong$bmi[!is.na(deathdatasetlong$bmi)], 50)

table(deathdatasetlong$bmi_cat, useNA = "always")
round(prop.table(table(deathdatasetlong$bmi_cat)) * 100, 1)

# ------------------------------------------------------------------------------
# 5. CREATE TWO COHORTS: 1996-2006 and 2008-2018
# ------------------------------------------------------------------------------

# Cohort 1: 1996-2006
cohort_1996 <- deathdatasetlong %>%
  filter(wave > 1994 & wave < 2008) %>%
  mutate(
    survival_timeYear = radyear - wave,
    survival_timeYear = ifelse(radyear > 2006 | is.na(radyear), 11, survival_timeYear),
    survival_timeMonth = ifelse(radyear > 2006 | is.na(radyear), 121, Survivalmonth),
    survival_timeMonth = pmin(survival_timeMonth, 121),
    # Censor at end of follow-up
    event = ifelse(survival_timeMonth == 121 & event == 1, 0, event),
    periodsince = "1996-2006"
  )

# Cohort 2: 2008-2018
cohort_2008 <- deathdatasetlong %>%
  filter(wave > 2006 & wave < 2020) %>%
  mutate(
    survival_timeYear = radyear - wave,
    survival_timeYear = ifelse(radyear > 2018 | is.na(radyear), 11, survival_timeYear),
    survival_timeMonth = ifelse(radyear > 2018 | is.na(radyear), 121, Survivalmonth),
    survival_timeMonth = pmin(survival_timeMonth, 121),
    event = ifelse(survival_timeMonth == 121 & event == 1, 0, event),
    periodsince = "2008-2018"
  ) %>%
  filter(survival_timeYear >= 0 | is.na(survival_timeYear))

# Combine cohorts
twosetofHRS <- rbind(cohort_1996, cohort_2008)

# ------------------------------------------------------------------------------
# 6. CREATE EDUCATION CATEGORIES
# ------------------------------------------------------------------------------

twosetofHRS <- twosetofHRS %>%
  mutate(
    education_level = case_when(
      raeduc %in% c(1, 2) ~ "Low",      # HS or less
      raeduc %in% c(3, 4) ~ "Medium",   # HS grad, some college
      raeduc == 5 ~ "High",             # College+
      TRUE ~ NA_character_
    ),
    education_level = factor(education_level, levels = c("Low", "Medium", "High"))
  )

# ------------------------------------------------------------------------------
# 7. CREATE STROKE VARIABLES (BASELINE, INCIDENT, RECURRENT)
# ------------------------------------------------------------------------------

twosetofHRSedit <- twosetofHRS %>%
  arrange(hhidpn, wave) %>%
  group_by(hhidpn) %>%
  mutate(
    # Identify cohort membership
    cohort_1996 = any(wave <= 2006),
    cohort_2008 = any(wave >= 2008),
    
    # Baseline waves
    baseline_wave_1996 = ifelse(cohort_1996, min(wave[wave <= 2006]), NA),
    baseline_wave_2008 = ifelse(cohort_2008, min(wave[wave >= 2008]), NA),
    
    # Is this the baseline wave?
    is_baseline_1996 = (wave == baseline_wave_1996 & !is.na(baseline_wave_1996)),
    is_baseline_2008 = (wave == baseline_wave_2008 & !is.na(baseline_wave_2008)),
    
    # Baseline stroke status for each cohort
    stroke_baseline_1996 = case_when(
      is_baseline_1996 & strok == "Stroke" ~ 1,
      is_baseline_1996 & strok == "No Stroke" ~ 0,
      TRUE ~ NA_real_
    ),
    
    stroke_baseline_2008 = case_when(
      is_baseline_2008 & strok == "Stroke" ~ 1,
      is_baseline_2008 & strok == "No Stroke" ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Fill forward within cohort
    stroke_baseline_1996 = ifelse(cohort_1996, max(stroke_baseline_1996, na.rm = TRUE), NA_real_),
    stroke_baseline_2008 = ifelse(cohort_2008, max(stroke_baseline_2008, na.rm = TRUE), NA_real_),
    
    stroke_baseline_1996 = ifelse(is.infinite(stroke_baseline_1996), NA_real_, stroke_baseline_1996),
    stroke_baseline_2008 = ifelse(is.infinite(stroke_baseline_2008), NA_real_, stroke_baseline_2008),
    
    # Combined baseline
    stroke_baseline = case_when(
      wave <= 2006 ~ stroke_baseline_1996,
      wave >= 2008 ~ stroke_baseline_2008,
      TRUE ~ NA_real_
    ),
    
    # Incident stroke (first stroke during follow-up)
    stroke_incident = case_when(
      is.na(stroke_baseline) | is.na(stroks) ~ NA_real_,
      stroke_baseline == 0 & stroks == 1 ~ 1,
      TRUE ~ 0
    ),
    
    # Recurrent stroke (stroke in those with baseline history)
    stroke_recurrent = case_when(
      is.na(stroke_baseline) | is.na(stroks) ~ NA_real_,
      stroke_baseline == 1 & stroks == 1 ~ 1,
      TRUE ~ 0
    ),
    
    # Combined acute stroke (incident OR recurrent)
    stroke_acutecombo = case_when(
      is.na(stroke_incident) | is.na(stroke_recurrent) ~ NA_real_,
      stroke_incident == 1 | stroke_recurrent == 1 ~ 1,
      TRUE ~ 0
    )
  ) %>%
  ungroup()

# ------------------------------------------------------------------------------
# 8. VERIFY DATA QUALITY
# ------------------------------------------------------------------------------

cat("\n=== DATA VERIFICATION ===\n")
cat("Total observations:", nrow(twosetofHRSedit), "\n")
cat("Unique individuals:", n_distinct(twosetofHRSedit$hhidpn), "\n\n")

cat("Baseline stroke by period:\n")
print(table(twosetofHRSedit$stroke_baseline, twosetofHRSedit$periodsince, useNA = "always"))

cat("\nAcute stroke by period:\n")
print(table(twosetofHRSedit$stroke_acutecombo, twosetofHRSedit$periodsince, useNA = "always"))

cat("\nEvents by period:\n")
print(table(twosetofHRSedit$event, twosetofHRSedit$periodsince, useNA = "always"))

# Check other variables
cat("\nRace/ethnicity by period:\n")
print(table(twosetofHRSedit$race_eth, twosetofHRSedit$periodsince, useNA = "always"))

cat("\nBMI category by period:\n")
print(table(twosetofHRSedit$bmi_cat, twosetofHRSedit$periodsince, useNA = "always"))

cat("\nInsurance type by period:\n")
print(table(twosetofHRSedit$insurance_type, twosetofHRSedit$periodsince, useNA = "always"))

cat("\nHeart disease by period:\n")
print(table(twosetofHRSedit$heart, twosetofHRSedit$periodsince, useNA = "always"))

# Check missingness for all model variables
model_vars <- c("survival_timeYear", "event", "stroke_baseline",
                "stroke_acutecombo", "periodsince", "ragender",
                "Age", "education_level", "smokev", "hibpe", "diabe",
                "heart", "bmi_cat", "insurance_type", "race_eth")
cat("\nMissing data by variable:\n")
print(colSums(is.na(twosetofHRSedit[, model_vars])))

# ==============================================================================
# 9. TABLE 1 
# ==============================================================================


# Create Table 1 variables
table1_vars <- c("Age", "education_level", "smokev", "hibpe", "diabe",
                 "heart", "bmi_cat", "insurance_type", "race_eth",
                 "stroke_baseline", "stroke_acutecombo", "event")

# Use baseline observations only for Table 1 (first observation per person per period)
table1_data <- twosetofHRSedit %>%
  group_by(hhidpn, periodsince) %>%
  filter(wave == min(wave)) %>%
  ungroup() %>%
  mutate(
    smokev_f = factor(smokev, levels = c(0, 1), labels = c("No", "Yes")),
    hibpe_f = factor(hibpe, levels = c(0, 1), labels = c("No", "Yes")),
    diabe_f = factor(diabe, levels = c(0, 1), labels = c("No", "Yes")),
    heart_f = factor(heart, levels = c(0, 1), labels = c("No", "Yes")),
    stroke_baseline_f = factor(stroke_baseline, levels = c(0, 1),
                               labels = c("No", "Yes")),
    event_f = factor(event, levels = c(0, 1), labels = c("No", "Yes"))
  )

# Table 1 stratified by period
cat("\n--- TABLE 1: By Period ---\n")
table1_by_period <- CreateTableOne(
  vars = c("Age", "education_level", "smokev_f", "hibpe_f", "diabe_f",
           "heart_f", "bmi_cat", "insurance_type", "race_eth",
           "stroke_baseline_f"),
  strata = "periodsince",
  data = table1_data,
  test = TRUE
)
print(table1_by_period, showAllLevels = TRUE, smd = FALSE)

# Table 1 stratified by period AND sex
cat("\n--- TABLE 1: By Period and Sex ---\n")
table1_data <- table1_data %>%
  mutate(period_sex = paste(periodsince, ragender, sep = " - "))

table1_by_period_sex <- CreateTableOne(
  vars = c("Age", "education_level", "smokev_f", "hibpe_f", "diabe_f",
           "heart_f", "bmi_cat", "insurance_type", "race_eth",
           "stroke_baseline_f"),
  strata = "period_sex",
  data = table1_data,
  test = TRUE
)
print(table1_by_period_sex, showAllLevels = TRUE, smd = FALSE)

# Incidence and prevalence rates per person-year
cat("\n--- STROKE RATES PER PERSON-YEAR ---\n")

rates_by_period <- twosetofHRSedit %>%
  group_by(periodsince, ragender) %>%
  summarise(
    n_individuals = n_distinct(hhidpn),
    total_person_years = sum(survival_timeYear, na.rm = TRUE),
    n_deaths = sum(event, na.rm = TRUE),
    n_baseline_stroke = sum(stroke_baseline == 1, na.rm = TRUE),
    n_incident_stroke = sum(stroke_acutecombo == 1, na.rm = TRUE),
    mortality_rate_per_1000py = round(n_deaths / total_person_years * 1000, 2),
    stroke_incidence_per_1000py = round(n_incident_stroke / total_person_years * 1000, 2),
    .groups = "drop"
  )

print(rates_by_period)

# Overall rates by period (both sexes)
rates_by_period_overall <- twosetofHRSedit %>%
  group_by(periodsince) %>%
  summarise(
    n_individuals = n_distinct(hhidpn),
    total_person_years = sum(survival_timeYear, na.rm = TRUE),
    n_deaths = sum(event, na.rm = TRUE),
    n_baseline_stroke = sum(stroke_baseline == 1, na.rm = TRUE),
    n_incident_stroke = sum(stroke_acutecombo == 1, na.rm = TRUE),
    mortality_rate_per_1000py = round(n_deaths / total_person_years * 1000, 2),
    stroke_incidence_per_1000py = round(n_incident_stroke / total_person_years * 1000, 2),
    .groups = "drop"
  )

cat("\nOverall rates by period:\n")
print(rates_by_period_overall)



# BMI by period and sex (4 columns for Table 1)
table(table1_data$bmi_cat, paste(table1_data$periodsince, table1_data$ragender))

# With percentages within each column
round(prop.table(table(table1_data$bmi_cat, 
                       paste(table1_data$periodsince, table1_data$ragender)), 
                 margin = 2) * 100, 1)


# Continuous BMI mean and SD by period and sex
table1_data %>%
  group_by(periodsince, ragender) %>%
  summarise(
    n = sum(!is.na(bmi_real)),
    mean_bmi = round(mean(bmi_real, na.rm = TRUE), 1),
    sd_bmi = round(sd(bmi_real, na.rm = TRUE), 1),
    .groups = "drop"
  )


table1_data %>%
  group_by(periodsince, ragender) %>%
  summarise(
    n = n(),
    age_mean = round(mean(Age, na.rm = TRUE), 1),
    age_sd = round(sd(Age, na.rm = TRUE), 1),
    stroke_mean = round(mean(stroke_baseline, na.rm = TRUE), 3),
    stroke_sd = round(sd(stroke_baseline, na.rm = TRUE), 3),
    smokev_mean = round(mean(smokev, na.rm = TRUE), 3),
    smokev_sd = round(sd(smokev, na.rm = TRUE), 3),
    hibpe_mean = round(mean(hibpe, na.rm = TRUE), 3),
    hibpe_sd = round(sd(hibpe, na.rm = TRUE), 3),
    diabe_mean = round(mean(diabe, na.rm = TRUE), 3),
    diabe_sd = round(sd(diabe, na.rm = TRUE), 3),
    heart_mean = round(mean(heart, na.rm = TRUE), 3),
    heart_sd = round(sd(heart, na.rm = TRUE), 3),
    bmi_mean = round(mean(bmi_real, na.rm = TRUE), 1),
    bmi_sd = round(sd(bmi_real, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  print(width = Inf)

# ==============================================================================
# 10. COX PROPORTIONAL HAZARDS MODELS - STRATIFIED BY SEX
# ==============================================================================

# ==============================================================================
# MALES
# ==============================================================================

cat("\n========== MALES ==========\n")

cat("\n=== MODEL 1: MALES ===\n")
cox1_male <- coxph(Surv(survival_timeYear, event) ~
                     stroke_baseline + stroke_baseline:periodsince +
                     periodsince + Age,
                   data = filter(twosetofHRSedit, ragender == "Male"))
print(summary(cox1_male))
cat("PH Test:\n")
print(cox.zph(cox1_male))

cat("\n=== MODEL 2: MALES ===\n")
cox2_male <- coxph(Surv(survival_timeYear, event) ~
                     stroke_baseline + stroke_acutecombo +
                     stroke_baseline:periodsince + stroke_acutecombo:periodsince +
                     periodsince + Age,
                   data = filter(twosetofHRSedit, ragender == "Male"))
print(summary(cox2_male))
cat("PH Test:\n")
print(cox.zph(cox2_male))

# cat("\n=== MODEL 3: MALES ===\n")
# cox3_male <- coxph(Surv(survival_timeYear, event) ~
#                      stroke_baseline + stroke_acutecombo +
#                      stroke_baseline:periodsince + stroke_acutecombo:periodsince +
#                      periodsince + Age + education_level + smokev + hibpe + diabe,
#                    data = filter(twosetofHRSedit, ragender == "Male"))
# print(summary(cox3_male))
# cat("PH Test:\n")
# print(cox.zph(cox3_male))
# 
# cat("\n=== MODEL 4: MALES ===\n")
# cox4_male <- coxph(Surv(survival_timeYear, event) ~
#                      stroke_baseline + stroke_acutecombo +
#                      stroke_baseline:periodsince + stroke_acutecombo:periodsince +
#                      periodsince + Age + education_level + smokev + hibpe + diabe +
#                      smokev:periodsince + education_level:periodsince +
#                      hibpe:periodsince + diabe:periodsince,
#                    data = filter(twosetofHRSedit, ragender == "Male"))
# print(summary(cox4_male))
# cat("PH Test:\n")
# print(cox.zph(cox4_male))

#  MODEL 3 - EXPANDED (add heart disease, BMI, insurance, race/ethnicity)
cat("\n=== MODEL 3: MALES (EXPANDED - additional confounders) ===\n")
cox3_male <- coxph(Surv(survival_timeYear, event) ~
                     stroke_baseline + stroke_acutecombo +
                     stroke_baseline:periodsince + stroke_acutecombo:periodsince +
                     periodsince + Age + education_level + smokev + hibpe + diabe +
                     heart + bmi_cat + insurance_type + race_eth,
                   data = filter(twosetofHRSedit, ragender == "Male"))
print(summary(cox3_male))
cat("PH Test:\n")
print(cox.zph(cox3_male))

#  MODEL 4 - EXPANDED + PERIOD INTERACTIONS for new covariates
cat("\n=== MODEL 4: MALES (EXPANDED + period interactions) ===\n")
cox4_male <- coxph(Surv(survival_timeYear, event) ~
                     stroke_baseline + stroke_acutecombo +
                     stroke_baseline:periodsince + stroke_acutecombo:periodsince +
                     periodsince + Age + education_level + smokev + hibpe + diabe +
                     heart + bmi_cat + insurance_type + race_eth +
                     smokev:periodsince + education_level:periodsince +
                     hibpe:periodsince + diabe:periodsince +
                     heart:periodsince + bmi_cat:periodsince +
                     insurance_type:periodsince + race_eth:periodsince,
                   data = filter(twosetofHRSedit, ragender == "Male"))
print(summary(cox4_male))
cat("PH Test:\n")
print(cox.zph(cox4_male))

# ==============================================================================
# FEMALES
# ==============================================================================

cat("\n========== FEMALES ==========\n")

cat("\n=== MODEL 1: FEMALES ===\n")
cox1_female <- coxph(Surv(survival_timeYear, event) ~
                       stroke_baseline + stroke_baseline:periodsince +
                       periodsince + Age,
                     data = filter(twosetofHRSedit, ragender == "Female"))
print(summary(cox1_female))
cat("PH Test:\n")
print(cox.zph(cox1_female))

cat("\n=== MODEL 2: FEMALES ===\n")
cox2_female <- coxph(Surv(survival_timeYear, event) ~
                       stroke_baseline + stroke_acutecombo +
                       stroke_baseline:periodsince + stroke_acutecombo:periodsince +
                       periodsince + Age,
                     data = filter(twosetofHRSedit, ragender == "Female"))
print(summary(cox2_female))
cat("PH Test:\n")
print(cox.zph(cox2_female))
# 
# cat("\n=== MODEL 3: FEMALES ===\n")
# cox3_female <- coxph(Surv(survival_timeYear, event) ~
#                        stroke_baseline + stroke_acutecombo +
#                        stroke_baseline:periodsince + stroke_acutecombo:periodsince +
#                        periodsince + Age + education_level + smokev + hibpe + diabe,
#                      data = filter(twosetofHRSedit, ragender == "Female"))
# print(summary(cox3_female))
# cat("PH Test:\n")
# print(cox.zph(cox3_female))
# 
# cat("\n=== MODEL 4: FEMALES ===\n")
# cox4_female <- coxph(Surv(survival_timeYear, event) ~
#                        stroke_baseline + stroke_acutecombo +
#                        stroke_baseline:periodsince + stroke_acutecombo:periodsince +
#                        periodsince + Age + education_level + smokev + hibpe + diabe +
#                        smokev:periodsince + education_level:periodsince +
#                        hibpe:periodsince + diabe:periodsince,
#                      data = filter(twosetofHRSedit, ragender == "Female"))
# print(summary(cox4_female))
# cat("PH Test:\n")
# print(cox.zph(cox4_female))

#  MODEL 3 - EXPANDED
cat("\n=== MODEL 3: FEMALES (EXPANDED - additional confounders) ===\n")
cox3_female <- coxph(Surv(survival_timeYear, event) ~
                       stroke_baseline + stroke_acutecombo +
                       stroke_baseline:periodsince + stroke_acutecombo:periodsince +
                       periodsince + Age + education_level + smokev + hibpe + diabe +
                       heart + bmi_cat + insurance_type + race_eth,
                     data = filter(twosetofHRSedit, ragender == "Female"))
print(summary(cox3_female))
cat("PH Test:\n")
print(cox.zph(cox3_female))

# MODEL 4 - EXPANDED + PERIOD INTERACTIONS
cat("\n=== MODEL 4: FEMALES (EXPANDED + period interactions) ===\n")
cox4_female <- coxph(Surv(survival_timeYear, event) ~
                       stroke_baseline + stroke_acutecombo +
                       stroke_baseline:periodsince + stroke_acutecombo:periodsince +
                       periodsince + Age + education_level + smokev + hibpe + diabe +
                       heart + bmi_cat + insurance_type + race_eth +
                       smokev:periodsince + education_level:periodsince +
                       hibpe:periodsince + diabe:periodsince +
                       heart:periodsince + bmi_cat:periodsince +
                       insurance_type:periodsince + race_eth:periodsince,
                     data = filter(twosetofHRSedit, ragender == "Female"))
print(summary(cox4_female))
cat("PH Test:\n")
print(cox.zph(cox4_female))

# ==============================================================================
# SCHOENFELD RESIDUAL TESTS SUMMARY
# ==============================================================================

cat("\n=== SCHOENFELD TESTS - ALL MODELS ===\n")
cat("\nMales:\n")
print(cox.zph(cox1_male))
print(cox.zph(cox2_male))
print(cox.zph(cox3_male))
print(cox.zph(cox4_male))

cat("\nFemales:\n")
print(cox.zph(cox1_female))
print(cox.zph(cox2_female))
print(cox.zph(cox3_female))
print(cox.zph(cox4_female))

# ==============================================================================
# MISSING DATA REPORT
# ==============================================================================

cat("\n=== MISSING DATA REPORT ===\n")

# Males
male_data <- filter(twosetofHRSedit, ragender == "Male")
total_males <- male_data %>% distinct(hhidpn) %>% nrow()

complete_males_m4 <- male_data %>%
  filter(!is.na(stroke_baseline) & !is.na(stroke_acutecombo) &
           !is.na(Age) & !is.na(education_level) &
           !is.na(smokev) & !is.na(hibpe) & !is.na(diabe) &
           !is.na(survival_timeYear) & !is.na(event)) %>%
  distinct(hhidpn) %>% nrow()

complete_males_m5 <- male_data %>%
  filter(!is.na(stroke_baseline) & !is.na(stroke_acutecombo) &
           !is.na(Age) & !is.na(education_level) &
           !is.na(smokev) & !is.na(hibpe) & !is.na(diabe) &
           !is.na(heart) & !is.na(bmi_cat) & !is.na(insurance_type) &
           !is.na(race_eth) &
           !is.na(survival_timeYear) & !is.na(event)) %>%
  distinct(hhidpn) %>% nrow()

cat("Males:\n")
cat("Total individuals:", total_males, "\n")
cat("In Models 1-4:", complete_males_m4, "\n")

cat("\nMales - Individuals with missing values:\n")
cat("Education:", male_data %>% filter(is.na(education_level)) %>% distinct(hhidpn) %>% nrow(), "\n")
cat("Smoking:", male_data %>% filter(is.na(smokev)) %>% distinct(hhidpn) %>% nrow(), "\n")
cat("Hypertension:", male_data %>% filter(is.na(hibpe)) %>% distinct(hhidpn) %>% nrow(), "\n")
cat("Diabetes:", male_data %>% filter(is.na(diabe)) %>% distinct(hhidpn) %>% nrow(), "\n")
cat("Heart disease:", male_data %>% filter(is.na(heart)) %>% distinct(hhidpn) %>% nrow(), "\n")
cat("BMI:", male_data %>% filter(is.na(bmi_cat)) %>% distinct(hhidpn) %>% nrow(), "\n")
cat("Insurance:", male_data %>% filter(is.na(insurance_type)) %>% distinct(hhidpn) %>% nrow(), "\n")
cat("Race/ethnicity:", male_data %>% filter(is.na(race_eth)) %>% distinct(hhidpn) %>% nrow(), "\n")

# Females
female_data <- filter(twosetofHRSedit, ragender == "Female")
total_females <- female_data %>% distinct(hhidpn) %>% nrow()

complete_females_m4 <- female_data %>%
  filter(!is.na(stroke_baseline) & !is.na(stroke_acutecombo) &
           !is.na(Age) & !is.na(education_level) &
           !is.na(smokev) & !is.na(hibpe) & !is.na(diabe) &
           !is.na(survival_timeYear) & !is.na(event)) %>%
  distinct(hhidpn) %>% nrow()

complete_females_m5 <- female_data %>%
  filter(!is.na(stroke_baseline) & !is.na(stroke_acutecombo) &
           !is.na(Age) & !is.na(education_level) &
           !is.na(smokev) & !is.na(hibpe) & !is.na(diabe) &
           !is.na(heart) & !is.na(bmi_cat) & !is.na(insurance_type) &
           !is.na(race_eth) &
           !is.na(survival_timeYear) & !is.na(event)) %>%
  distinct(hhidpn) %>% nrow()

cat("\nFemales:\n")
cat("Total individuals:", total_females, "\n")
cat("In Models 1-4:", complete_females_m4, "\n")

cat("\nFemales - Individuals with missing values:\n")
cat("Education:", female_data %>% filter(is.na(education_level)) %>% distinct(hhidpn) %>% nrow(), "\n")
cat("Smoking:", female_data %>% filter(is.na(smokev)) %>% distinct(hhidpn) %>% nrow(), "\n")
cat("Hypertension:", female_data %>% filter(is.na(hibpe)) %>% distinct(hhidpn) %>% nrow(), "\n")
cat("Diabetes:", female_data %>% filter(is.na(diabe)) %>% distinct(hhidpn) %>% nrow(), "\n")
cat("Heart disease:", female_data %>% filter(is.na(heart)) %>% distinct(hhidpn) %>% nrow(), "\n")
cat("BMI:", female_data %>% filter(is.na(bmi_cat)) %>% distinct(hhidpn) %>% nrow(), "\n")
cat("Insurance:", female_data %>% filter(is.na(insurance_type)) %>% distinct(hhidpn) %>% nrow(), "\n")
cat("Race/ethnicity:", female_data %>% filter(is.na(race_eth)) %>% distinct(hhidpn) %>% nrow(), "\n")

# Person-years
personyears_males <- male_data %>%
  filter(!is.na(stroke_baseline) & !is.na(stroke_acutecombo) &
           !is.na(Age) & !is.na(education_level) &
           !is.na(smokev) & !is.na(hibpe) & !is.na(diabe) &
           !is.na(survival_timeYear) & !is.na(event)) %>%
  nrow()

personyears_females <- female_data %>%
  filter(!is.na(stroke_baseline) & !is.na(stroke_acutecombo) &
           !is.na(Age) & !is.na(education_level) &
           !is.na(smokev) & !is.na(hibpe) & !is.na(diabe) &
           !is.na(survival_timeYear) & !is.na(event)) %>%
  nrow()

cat("\nMales person-years (obs):", personyears_males, "\n")
cat("Females person-years (obs):", personyears_females, "\n")
cat("Total person-years (obs):", personyears_males + personyears_females, "\n")

# ==============================================================================
# 11. MODEL COMPARISON (AIC)
# ==============================================================================

cat("\n=== MODEL COMPARISON (AIC) - MALES ===\n")
cat("Model 1 AIC:", AIC(cox1_male), "\n")
cat("Model 2 AIC:", AIC(cox2_male), "\n")
cat("Model 3 AIC:", AIC(cox3_male), "\n")
cat("Model 4 AIC:", AIC(cox4_male), "\n")

cat("\n=== MODEL COMPARISON (AIC) - FEMALES ===\n")
cat("Model 1 AIC:", AIC(cox1_female), "\n")
cat("Model 2 AIC:", AIC(cox2_female), "\n")
cat("Model 3 AIC:", AIC(cox3_female), "\n")
cat("Model 4 AIC:", AIC(cox4_female), "\n")

# # Likelihood ratio tests for new models
# cat("\n=== LIKELIHOOD RATIO TESTS ===\n")
# cat("Males - Model 5 vs Model 3:\n")
# print(anova(cox3_male, cox3_male))
# cat("\nMales - Model 6 vs Model 5:\n")
# print(anova(cox3_male, cox4_male))
# 
# cat("\nFemales - Model 5 vs Model 3:\n")
# print(anova(cox3_female, cox3_female))
# cat("\nFemales - Model 6 vs Model 5:\n")
# print(anova(cox3_female, cox4_female))

# ==============================================================================
# 12. MODEL FIT STATISTICS
# ==============================================================================

calculate_rsquared <- function(cox_model) {
  loglik_null <- cox_model$loglik[1]
  loglik_model <- cox_model$loglik[2]
  n <- cox_model$n
  
  nagelkerke_r2 <- (1 - exp((2/n) * (loglik_null - loglik_model))) /
    (1 - exp((2/n) * loglik_null))
  cox_snell_r2 <- 1 - exp((2/n) * (loglik_null - loglik_model))
  
  return(list(
    nagelkerke = nagelkerke_r2,
    cox_snell = cox_snell_r2,
    concordance = cox_model$concordance['concordance'],
    n = n,
    events = cox_model$nevent,
    aic = AIC(cox_model)
  ))
}

cat("\n=== MODEL FIT STATISTICS - MALES ===\n")
for(i in 1:4) {
  model <- get(paste0("cox", i, "_male"))
  r2 <- calculate_rsquared(model)
  cat(sprintf("Model %d: N=%d, Events=%d, AIC=%.1f, C=%.3f, Nagelkerke R2=%.3f\n",
              i, r2$n, r2$events, r2$aic, r2$concordance, r2$nagelkerke))
}

cat("\n=== MODEL FIT STATISTICS - FEMALES ===\n")
for(i in 1:4) {
  model <- get(paste0("cox", i, "_female"))
  r2 <- calculate_rsquared(model)
  cat(sprintf("Model %d: N=%d, Events=%d, AIC=%.1f, C=%.3f, Nagelkerke R2=%.3f\n",
              i, r2$n, r2$events, r2$aic, r2$concordance, r2$nagelkerke))
}

# ==============================================================================
# 13. DTMS ANALYSIS MODELS 
# ==============================================================================

cat("\n=== PREPARING DATA FOR DTMS ===\n")

# Create 1996-2006 cohort for DTMS
onlyHRSsince1996 <- twosetofHRSedit %>%
  filter(periodsince == "1996-2006") %>%
  filter(Age > 49.99 & Age < 109) %>%
  arrange(hhidpn, wave) %>%
  group_by(hhidpn) %>%
  mutate(
    Age2 = round(Age, 0),
    ever_had_stroke = cummax(ifelse(
      stroke_baseline == 1 | (stroke_acutecombo == 1 & !is.na(stroke_acutecombo)),
      1, 0
    )),
    state_dtms = case_when(
      event == 1 ~ "Dead",
      stroke_baseline == 1 | stroke_acutecombo == 1 | ever_had_stroke == 1 ~ "Stroke History",
      ever_had_stroke == 0 ~ "No Stroke",
      is.na(stroke_baseline) & is.na(stroke_acutecombo) ~ "No Stroke",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

# Create 2008-2018 cohort for DTMS
onlyHRSsince2008 <- twosetofHRSedit %>%
  filter(periodsince == "2008-2018") %>%
  filter(Age > 49.99 & Age < 109) %>%
  arrange(hhidpn, wave) %>%
  group_by(hhidpn) %>%
  mutate(
    Age2 = round(Age, 0),
    ever_had_stroke = cummax(ifelse(
      stroke_baseline == 1 | (stroke_acutecombo == 1 & !is.na(stroke_acutecombo)),
      1, 0
    )),
    state_dtms = case_when(
      event == 1 ~ "Dead",
      stroke_baseline == 1 | stroke_acutecombo == 1 | ever_had_stroke == 1 ~ "Stroke History",
      ever_had_stroke == 0 ~ "No Stroke",
      is.na(stroke_baseline) & is.na(stroke_acutecombo) ~ "No Stroke",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

# Verify state distributions
cat("\nState distribution 1996-2006:\n")
print(table(onlyHRSsince1996$state_dtms, useNA = "always"))
cat("\nState distribution 2008-2018:\n")
print(table(onlyHRSsince2008$state_dtms, useNA = "always"))

# DTMS model definition
definitiondtms1996 <- dtms(
  transient = c("No Stroke", "Stroke History"),
  absorbing = "Dead",
  timescale = 50:108,
  timestep = 2
)

# 1996-2006 FULL SAMPLE
Modelfor1996 <- dtms_format(
  data = onlyHRSsince1996, dtms = definitiondtms1996,
  idvar = "hhidpn", timevar = "Age2", statevar = "state_dtms"
)
Modelfor1996 <- dtms_clean(data = Modelfor1996, dtms = definitiondtms1996)
fit1996 <- dtms_fit(data = Modelfor1996, controls = "time")
probs1996 <- dtms_transitions(model = fit1996, dtms = definitiondtms1996)
Tp1996 <- dtms_matrix(probs = probs1996, dtms = definitiondtms1996)
S1996 <- dtms_start(dtms = definitiondtms1996, data = Modelfor1996, start_time = 50)

LE1996 <- dtms_expectancy(
  dtms = definitiondtms1996, matrix = Tp1996,
  start_distr = S1996, start_time = 50
)
cat("\nLife expectancies at age 50 (1996-2006):\n")
print(LE1996)

bootfun <- function(data, dtms) {
  fit <- dtms_fit(data = data, controls = "time")
  probs <- dtms_transitions(dtms = dtms, model = fit)
  Tp <- dtms_matrix(dtms = dtms, probs = probs)
  S <- dtms_start(dtms = dtms, data = data)
  dtms_expectancy(dtms = dtms, matrix = Tp, start_distr = S)
}

cat("\nBootstrapping 1996-2006 (50 replications)...\n")
bootresults1996 <- dtms_boot(
  data = Modelfor1996, dtms = definitiondtms1996,
  fun = bootfun, idvar = "id", rep = 1000,
  method = "block", parallel = TRUE
)
cat("\n1996-2006 Bootstrap Results:\n")
print(summary(bootresults1996))

# 2008-2018 FULL SAMPLE
definitiondtms2008 <- dtms(
  transient = c("No Stroke", "Stroke History"),
  absorbing = "Dead",
  timescale = 50:108,
  timestep = 2
)

Modelfor2008 <- dtms_format(
  data = onlyHRSsince2008, dtms = definitiondtms2008,
  idvar = "hhidpn", timevar = "Age2", statevar = "state_dtms"
)
Modelfor2008 <- dtms_clean(data = Modelfor2008, dtms = definitiondtms2008)
fit2008 <- dtms_fit(data = Modelfor2008, controls = "time")
probs2008 <- dtms_transitions(model = fit2008, dtms = definitiondtms2008)
Tp2008 <- dtms_matrix(probs = probs2008, dtms = definitiondtms2008)
S2008 <- dtms_start(dtms = definitiondtms2008, data = Modelfor2008, start_time = 50)

LE2008 <- dtms_expectancy(
  dtms = definitiondtms2008, matrix = Tp2008,
  start_distr = S2008, start_time = 50
)
cat("\nLife expectancies at age 50 (2008-2018):\n")
print(LE2008)

cat("\nBootstrapping 2008-2018 (50 replications)...\n")
bootresults2008 <- dtms_boot(
  data = Modelfor2008, dtms = definitiondtms2008,
  fun = bootfun, idvar = "id", rep = 1000,
  method = "block", parallel = TRUE
)
cat("\n2008-2018 Bootstrap Results:\n")
print(summary(bootresults2008))

# DTMS BY SEX - 1996-2006 MALES
cat("\n=== DTMS ANALYSIS: 1996-2006 MALES ===\n")
onlyHRSsince1996Male <- onlyHRSsince1996 %>% filter(ragender == "Male")
Modelfor1996Male <- dtms_format(
  data = onlyHRSsince1996Male, dtms = definitiondtms1996,
  idvar = "hhidpn", timevar = "Age2", statevar = "state_dtms"
)
Modelfor1996Male <- dtms_clean(data = Modelfor1996Male, dtms = definitiondtms1996)
fit1996Male <- dtms_fit(data = Modelfor1996Male, controls = "time")
probs1996Male <- dtms_transitions(model = fit1996Male, dtms = definitiondtms1996)
Tp1996Male <- dtms_matrix(probs = probs1996Male, dtms = definitiondtms1996)
S1996Male <- dtms_start(dtms = definitiondtms1996, data = Modelfor1996Male, start_time = 50)
LE1996Male <- dtms_expectancy(
  dtms = definitiondtms1996, matrix = Tp1996Male,
  start_distr = S1996Male, start_time = 50
)
cat("\nLife expectancies at age 50 (1996-2006 Males):\n")
print(LE1996Male)
cat("\nBootstrapping...\n")
bootresults1996Male <- dtms_boot(
  data = Modelfor1996Male, dtms = definitiondtms1996,
  fun = bootfun, idvar = "id", rep = 1000,
  method = "block", parallel = TRUE
)
cat("\n1996-2006 Males Bootstrap Results:\n")
print(summary(bootresults1996Male))

# DTMS BY SEX - 1996-2006 FEMALES
cat("\n=== DTMS ANALYSIS: 1996-2006 FEMALES ===\n")
onlyHRSsince1996Female <- onlyHRSsince1996 %>% filter(ragender == "Female")
Modelfor1996Female <- dtms_format(
  data = onlyHRSsince1996Female, dtms = definitiondtms1996,
  idvar = "hhidpn", timevar = "Age2", statevar = "state_dtms"
)
Modelfor1996Female <- dtms_clean(data = Modelfor1996Female, dtms = definitiondtms1996)
fit1996Female <- dtms_fit(data = Modelfor1996Female, controls = "time")
probs1996Female <- dtms_transitions(model = fit1996Female, dtms = definitiondtms1996)
Tp1996Female <- dtms_matrix(probs = probs1996Female, dtms = definitiondtms1996)
S1996Female <- dtms_start(dtms = definitiondtms1996, data = Modelfor1996Female, start_time = 50)
LE1996Female <- dtms_expectancy(
  dtms = definitiondtms1996, matrix = Tp1996Female,
  start_distr = S1996Female, start_time = 50
)
cat("\nLife expectancies at age 50 (1996-2006 Females):\n")
print(LE1996Female)
cat("\nBootstrapping...\n")
bootresults1996Female <- dtms_boot(
  data = Modelfor1996Female, dtms = definitiondtms1996,
  fun = bootfun, idvar = "id", rep = 1000,
  method = "block", parallel = TRUE
)
cat("\n1996-2006 Females Bootstrap Results:\n")
print(summary(bootresults1996Female))

# DTMS BY SEX - 2008-2018 MALES
cat("\n=== DTMS ANALYSIS: 2008-2018 MALES ===\n")
onlyHRSsince2008Male <- onlyHRSsince2008 %>% filter(ragender == "Male")
Modelfor2008Male <- dtms_format(
  data = onlyHRSsince2008Male, dtms = definitiondtms2008,
  idvar = "hhidpn", timevar = "Age2", statevar = "state_dtms"
)
Modelfor2008Male <- dtms_clean(data = Modelfor2008Male, dtms = definitiondtms2008)
fit2008Male <- dtms_fit(data = Modelfor2008Male, controls = "time")
probs2008Male <- dtms_transitions(model = fit2008Male, dtms = definitiondtms2008)
Tp2008Male <- dtms_matrix(probs = probs2008Male, dtms = definitiondtms2008)
S2008Male <- dtms_start(dtms = definitiondtms2008, data = Modelfor2008Male, start_time = 50)
LE2008Male <- dtms_expectancy(
  dtms = definitiondtms2008, matrix = Tp2008Male,
  start_distr = S2008Male, start_time = 50
)
cat("\nLife expectancies at age 50 (2008-2018 Males):\n")
print(LE2008Male)
cat("\nBootstrapping...\n")
bootresults2008Male <- dtms_boot(
  data = Modelfor2008Male, dtms = definitiondtms2008,
  fun = bootfun, idvar = "id", rep = 1000,
  method = "block", parallel = TRUE
)
cat("\n2008-2018 Males Bootstrap Results:\n")
print(summary(bootresults2008Male))

# DTMS BY SEX - 2008-2018 FEMALES
cat("\n=== DTMS ANALYSIS: 2008-2018 FEMALES ===\n")
onlyHRSsince2008Female <- onlyHRSsince2008 %>% filter(ragender == "Female")
Modelfor2008Female <- dtms_format(
  data = onlyHRSsince2008Female, dtms = definitiondtms2008,
  idvar = "hhidpn", timevar = "Age2", statevar = "state_dtms"
)
Modelfor2008Female <- dtms_clean(data = Modelfor2008Female, dtms = definitiondtms2008)
fit2008Female <- dtms_fit(data = Modelfor2008Female, controls = "time")
probs2008Female <- dtms_transitions(model = fit2008Female, dtms = definitiondtms2008)
Tp2008Female <- dtms_matrix(probs = probs2008Female, dtms = definitiondtms2008)
S2008Female <- dtms_start(dtms = definitiondtms2008, data = Modelfor2008Female, start_time = 50)
LE2008Female <- dtms_expectancy(
  dtms = definitiondtms2008, matrix = Tp2008Female,
  start_distr = S2008Female, start_time = 50
)
cat("\nLife expectancies at age 50 (2008-2018 Females):\n")
print(LE2008Female)
cat("\nBootstrapping...\n")
bootresults2008Female <- dtms_boot(
  data = Modelfor2008Female, dtms = definitiondtms2008,
  fun = bootfun, idvar = "id", rep = 1000,
  method = "block", parallel = TRUE
)
cat("\n2008-2018 Females Bootstrap Results:\n")
print(summary(bootresults2008Female))

# Life expectancy summary
cat("\n=== SUMMARY OF LIFE EXPECTANCIES AT AGE 50 ===\n")

# ##############################################################################
# ##############################################################################
# 14. SENSITIVITY ANALYSIS: EXCLUDING OVERLAPPING INDIVIDUALS
# ##############################################################################
# ##############################################################################

# ==============================================================================
# 14. SENSITIVITY ANALYSIS: RANDOM PERIOD ASSIGNMENT (Supplementary Table S8)
# ==============================================================================
# 17,334 individuals contribute observations to both periods. To check whether
# this overlap drives our results, we randomly assign each overlapping
# individual to ONE period only (with equal probability) and rerun Cox Model 3.
# Random assignment preserves both period samples while ensuring no individual
# is counted twice.
# ==============================================================================

# Identify overlapping individuals
ids_1996 <- twosetofHRSedit %>%
  filter(periodsince == "1996-2006") %>%
  distinct(hhidpn) %>% pull(hhidpn)

ids_2008 <- twosetofHRSedit %>%
  filter(periodsince == "2008-2018") %>%
  distinct(hhidpn) %>% pull(hhidpn)

overlap_ids <- intersect(ids_1996, ids_2008)

cat("\nIndividuals in 1996-2006 only:", length(setdiff(ids_1996, ids_2008)), "\n")
cat("Individuals in 2008-2018 only:", length(setdiff(ids_2008, ids_1996)), "\n")
cat("Individuals in BOTH periods:  ", length(overlap_ids), "\n")
cat("Total unique individuals:     ", length(union(ids_1996, ids_2008)), "\n\n")

# Random assignment, fixed seed for reproducibility
set.seed(12345)
overlap_assignment <- data.frame(
  hhidpn      = overlap_ids,
  keep_period = sample(c("1996-2006", "2008-2018"),
                       length(overlap_ids), replace = TRUE)
)
cat("Assigned to 1996-2006:", sum(overlap_assignment$keep_period == "1996-2006"), "\n")
cat("Assigned to 2008-2018:", sum(overlap_assignment$keep_period == "2008-2018"), "\n\n")

# Build dataset where overlapping individuals appear in only one period
twosetofHRSedit_random <- twosetofHRSedit %>%
  left_join(overlap_assignment, by = "hhidpn") %>%
  filter(is.na(keep_period) | periodsince == keep_period) %>%
  select(-keep_period)

# Verify no individual appears in both periods now
check <- twosetofHRSedit_random %>%
  distinct(hhidpn, periodsince) %>%
  group_by(hhidpn) %>% filter(n() > 1) %>% nrow()
cat("Individuals still in both periods (should be 0):", check, "\n\n")

cat("Sample sizes after random assignment:\n")
print(twosetofHRSedit_random %>%
        group_by(periodsince) %>%
        summarise(n_obs = n(), n_individuals = n_distinct(hhidpn), .groups = "drop"))

# --- Cox Model 3 (expanded): Males ---
cat("\n=== SENSITIVITY: MODEL 3 MALES (random assignment) ===\n")
cox3_male_rand <- coxph(Surv(survival_timeYear, event) ~
                          stroke_baseline + stroke_acutecombo +
                          stroke_baseline:periodsince + stroke_acutecombo:periodsince +
                          periodsince + Age + education_level + smokev + hibpe + diabe +
                          heart + bmi_cat + insurance_type + race_eth,
                        data = filter(twosetofHRSedit_random, ragender == "Male"))
print(summary(cox3_male_rand))

# --- Cox Model 3 (expanded): Females ---
cat("\n=== SENSITIVITY: MODEL 3 FEMALES (random assignment) ===\n")
cox3_female_rand <- coxph(Surv(survival_timeYear, event) ~
                            stroke_baseline + stroke_acutecombo +
                            stroke_baseline:periodsince + stroke_acutecombo:periodsince +
                            periodsince + Age + education_level + smokev + hibpe + diabe +
                            heart + bmi_cat + insurance_type + race_eth,
                          data = filter(twosetofHRSedit_random, ragender == "Female"))
print(summary(cox3_female_rand))


# ==============================================================================
# 15. SENSITIVITY ANALYSIS: DTMS OVER FULL 1996-2018 PERIOD
# ==============================================================================
# Estimates LE50 over the full 1996-2018 period without splitting into two
# 10-year windows. Reported for completeness; main analysis uses two-period.
# ==============================================================================

fullperiod <- deathdatasetlong %>%
  filter(wave >= 1996 & wave <= 2018) %>%
  filter(Age > 49.99 & Age < 109) %>%
  arrange(hhidpn, wave) %>%
  group_by(hhidpn) %>%
  mutate(
    Age2 = round(Age, 0),
    ever_had_stroke = cummax(ifelse(strok == "Stroke", 1, 0)),
    state_dtms = case_when(
      event == 1            ~ "Dead",
      ever_had_stroke == 1  ~ "Stroke History",
      ever_had_stroke == 0  ~ "No Stroke",
      TRUE                  ~ NA_character_
    )
  ) %>%
  ungroup()

cat("\n=== FULL PERIOD 1996-2018 ===\n")
cat("Total observations:", nrow(fullperiod), "\n")
cat("Unique individuals:", n_distinct(fullperiod$hhidpn), "\n")
cat("\nState distribution:\n"); print(table(fullperiod$state_dtms, useNA = "always"))
cat("\nBy sex:\n"); print(table(fullperiod$state_dtms, fullperiod$ragender, useNA = "always"))

definition_full <- dtms(
  transient = c("No Stroke", "Stroke History"),
  absorbing = "Dead",
  timescale = 50:108,
  timestep  = 2
)

# Both sexes
Model_full <- dtms_format(data = fullperiod, dtms = definition_full,
                          idvar = "hhidpn", timevar = "Age2", statevar = "state_dtms")
Model_full <- dtms_clean(data = Model_full, dtms = definition_full)
fit_full   <- dtms_fit(data = Model_full, controls = "time")
probs_full <- dtms_transitions(model = fit_full, dtms = definition_full)
Tp_full    <- dtms_matrix(probs = probs_full, dtms = definition_full)
S_full     <- dtms_start(dtms = definition_full, data = Model_full, start_time = 50)
LE_full    <- dtms_expectancy(dtms = definition_full, matrix = Tp_full,
                              start_distr = S_full, start_time = 50)
cat("\nLE50 (full period 1996-2018, both sexes):\n"); print(LE_full)

cat("\nBootstrapping 1996-2018 both sexes (1000 reps)...\n")
boot_full <- dtms_boot(data = Model_full, dtms = definition_full,
                       fun = bootfun, idvar = "id", rep = 1000,
                       method = "block", parallel = TRUE)
print(summary(boot_full))

# Males
fullperiod_male <- fullperiod %>% filter(ragender == "Male")
Model_full_male <- dtms_format(data = fullperiod_male, dtms = definition_full,
                               idvar = "hhidpn", timevar = "Age2", statevar = "state_dtms")
Model_full_male <- dtms_clean(data = Model_full_male, dtms = definition_full)
fit_full_male   <- dtms_fit(data = Model_full_male, controls = "time")
probs_full_male <- dtms_transitions(model = fit_full_male, dtms = definition_full)
Tp_full_male    <- dtms_matrix(probs = probs_full_male, dtms = definition_full)
S_full_male     <- dtms_start(dtms = definition_full, data = Model_full_male, start_time = 50)
LE_full_male    <- dtms_expectancy(dtms = definition_full, matrix = Tp_full_male,
                                   start_distr = S_full_male, start_time = 50)
cat("\nLE50 (full period, Males):\n"); print(LE_full_male)

cat("\nBootstrapping 1996-2018 Males (1000 reps)...\n")
boot_full_male <- dtms_boot(data = Model_full_male, dtms = definition_full,
                            fun = bootfun, idvar = "id", rep = 1000,
                            method = "block", parallel = TRUE)
print(summary(boot_full_male))

# Females
fullperiod_female <- fullperiod %>% filter(ragender == "Female")
Model_full_female <- dtms_format(data = fullperiod_female, dtms = definition_full,
                                 idvar = "hhidpn", timevar = "Age2", statevar = "state_dtms")
Model_full_female <- dtms_clean(data = Model_full_female, dtms = definition_full)
fit_full_female   <- dtms_fit(data = Model_full_female, controls = "time")
probs_full_female <- dtms_transitions(model = fit_full_female, dtms = definition_full)
Tp_full_female    <- dtms_matrix(probs = probs_full_female, dtms = definition_full)
S_full_female     <- dtms_start(dtms = definition_full, data = Model_full_female, start_time = 50)
LE_full_female    <- dtms_expectancy(dtms = definition_full, matrix = Tp_full_female,
                                     start_distr = S_full_female, start_time = 50)
cat("\nLE50 (full period, Females):\n"); print(LE_full_female)

cat("\nBootstrapping 1996-2018 Females (1000 reps)...\n")
boot_full_female <- dtms_boot(data = Model_full_female, dtms = definition_full,
                              fun = bootfun, idvar = "id", rep = 1000,
                              method = "block", parallel = TRUE)
print(summary(boot_full_female))


# ==============================================================================
# 16. SENSITIVITY ANALYSIS: COX MODELS WITH SEX INTERACTIONS (Table S9)
# ==============================================================================
# Tests whether (a) stroke mortality differs by sex and (b) the change in
# stroke mortality between periods differs by sex. Built on the Model 3
# specification (full set of confounders).
# ==============================================================================

# --- Stroke x Sex interaction (does stroke kill men and women equally?) ---
cat("\n=== Sex interactions: stroke x sex (Model 3 spec) ===\n")
cox_sex_main <- coxph(Surv(survival_timeYear, event) ~
                        stroke_baseline + stroke_acutecombo +
                        ragender + periodsince + Age +
                        education_level + smokev + hibpe + diabe +
                        heart + bmi_cat + insurance_type + race_eth +
                        stroke_baseline:periodsince + stroke_acutecombo:periodsince +
                        stroke_baseline:ragender + stroke_acutecombo:ragender,
                      data = twosetofHRSedit)
print(summary(cox_sex_main))
cat("\nSchoenfeld test:\n"); print(cox.zph(cox_sex_main))

cat("\nKey sex-related HRs:\n")
hr_main <- summary(cox_sex_main)$conf.int
print(round(hr_main[grep("ragender", rownames(hr_main)), ], 4))

# --- 3-way interaction (did men and women gain equally across periods?) ---
cat("\n=== Sex interactions: stroke x sex x period (3-way) ===\n")
cox_sex_3way <- coxph(Surv(survival_timeYear, event) ~
                        stroke_baseline + stroke_acutecombo +
                        ragender + periodsince + Age +
                        education_level + smokev + hibpe + diabe +
                        heart + bmi_cat + insurance_type + race_eth +
                        stroke_baseline:periodsince + stroke_acutecombo:periodsince +
                        stroke_baseline:ragender + stroke_acutecombo:ragender +
                        ragender:periodsince +
                        stroke_baseline:ragender:periodsince +
                        stroke_acutecombo:ragender:periodsince,
                      data = twosetofHRSedit)
print(summary(cox_sex_3way))
cat("\nSchoenfeld test:\n"); print(cox.zph(cox_sex_3way))

cat("\nKey sex-related HRs:\n")
hr_3way <- summary(cox_sex_3way)$conf.int
print(round(hr_3way[grep("ragender", rownames(hr_3way)), ], 4))

cat("\nModel comparison:\n")
cat("  AIC stroke x sex:           ", AIC(cox_sex_main), "\n")
cat("  AIC stroke x sex x period:  ", AIC(cox_sex_3way), "\n")
cat("\nLikelihood ratio test (3-way vs 2-way):\n")
print(anova(cox_sex_main, cox_sex_3way))

cat("\nINTERPRETATION:\n")
cat(" stroke_baseline:ragenderFemale\n")
cat("   HR < 1: women have LOWER stroke mortality than men.\n")
cat("   HR > 1: women have HIGHER stroke mortality than men.\n")
cat(" stroke_baseline:ragenderFemale:periodsince2008-2018\n")
cat("   Non-significant = gains were similar across sexes.\n")
cat("   Significant     = one sex gained more than the other.\n")


# ==============================================================================
# 17. SENSITIVITY ANALYSIS: DTMS WITH SEX AS COVARIATE (Table S11)
# ==============================================================================
# Three approaches to sex differences in DTMS:
#   (1) Separate models per sex/period (main analysis, Section 13).
#   (2) Sex as additive covariate: one DTMS per period, Gender shifts levels.
#   (3) Sex x Age interaction: one DTMS per period, Gender effect varies by age.
# All three approaches yield consistent sex differences. Approach 3 is the
# most realistic (female advantage shrinks at old ages).
# ==============================================================================

# Build datasets carrying numeric Gender into dtms_format
onlyHRSsince1996_withsex <- twosetofHRSedit %>%
  filter(periodsince == "1996-2006") %>%
  filter(Age > 49.99 & Age < 109) %>%
  arrange(hhidpn, wave) %>%
  group_by(hhidpn) %>%
  mutate(
    Age2   = round(Age, 0),
    Gender = ifelse(ragender == "Female", 1, 0),
    ever_had_stroke = cummax(ifelse(
      stroke_baseline == 1 | (stroke_acutecombo == 1 & !is.na(stroke_acutecombo)),
      1, 0)),
    state_dtms = case_when(
      event == 1                                        ~ "Dead",
      stroke_baseline == 1 | stroke_acutecombo == 1 |
        ever_had_stroke == 1                            ~ "Stroke History",
      ever_had_stroke == 0                              ~ "No Stroke",
      is.na(stroke_baseline) & is.na(stroke_acutecombo) ~ "No Stroke",
      TRUE                                              ~ NA_character_
    )
  ) %>% ungroup()

onlyHRSsince2008_withsex <- twosetofHRSedit %>%
  filter(periodsince == "2008-2018") %>%
  filter(Age > 49.99 & Age < 109) %>%
  arrange(hhidpn, wave) %>%
  group_by(hhidpn) %>%
  mutate(
    Age2   = round(Age, 0),
    Gender = ifelse(ragender == "Female", 1, 0),
    ever_had_stroke = cummax(ifelse(
      stroke_baseline == 1 | (stroke_acutecombo == 1 & !is.na(stroke_acutecombo)),
      1, 0)),
    state_dtms = case_when(
      event == 1                                        ~ "Dead",
      stroke_baseline == 1 | stroke_acutecombo == 1 |
        ever_had_stroke == 1                            ~ "Stroke History",
      ever_had_stroke == 0                              ~ "No Stroke",
      is.na(stroke_baseline) & is.na(stroke_acutecombo) ~ "No Stroke",
      TRUE                                              ~ NA_character_
    )
  ) %>% ungroup()

definitiondtms <- dtms(
  transient = c("No Stroke", "Stroke History"),
  absorbing = "Dead",
  timescale = 50:108,
  timestep  = 2
)

Model1996_sex <- dtms_format(data = onlyHRSsince1996_withsex, dtms = definitiondtms,
                             idvar = "hhidpn", timevar = "Age2", statevar = "state_dtms")
Model1996_sex <- dtms_clean(data = Model1996_sex, dtms = definitiondtms)

Model2008_sex <- dtms_format(data = onlyHRSsince2008_withsex, dtms = definitiondtms,
                             idvar = "hhidpn", timevar = "Age2", statevar = "state_dtms")
Model2008_sex <- dtms_clean(data = Model2008_sex, dtms = definitiondtms)

# --- Approach (2): Sex as additive covariate ---
cat("\n=== DTMS: sex as additive covariate ===\n")

fit1996_sex <- dtms_fit(data = Model1996_sex, controls = c("time", "Gender"))
fit2008_sex <- dtms_fit(data = Model2008_sex, controls = c("time", "Gender"))

# 1996-2006 Males
probs1996_male  <- dtms_transitions(model = fit1996_sex, dtms = definitiondtms,
                                    controls = list(time = definitiondtms$timescale, Gender = 0))
Tp1996_male     <- dtms_matrix(probs = probs1996_male, dtms = definitiondtms)
S1996_male      <- dtms_start(dtms = definitiondtms, data = Model1996_sex,
                              start_time = 50, variables = list(Gender = 0))
LE1996_male_sex <- dtms_expectancy(dtms = definitiondtms, matrix = Tp1996_male,
                                   start_distr = S1996_male, start_time = 50)
cat("\nLE50 1996-2006 Males (additive):\n"); print(LE1996_male_sex)

# 1996-2006 Females
probs1996_female  <- dtms_transitions(model = fit1996_sex, dtms = definitiondtms,
                                      controls = list(time = definitiondtms$timescale, Gender = 1))
Tp1996_female     <- dtms_matrix(probs = probs1996_female, dtms = definitiondtms)
S1996_female      <- dtms_start(dtms = definitiondtms, data = Model1996_sex,
                                start_time = 50, variables = list(Gender = 1))
LE1996_female_sex <- dtms_expectancy(dtms = definitiondtms, matrix = Tp1996_female,
                                     start_distr = S1996_female, start_time = 50)
cat("\nLE50 1996-2006 Females (additive):\n"); print(LE1996_female_sex)

# 2008-2018 Males
probs2008_male  <- dtms_transitions(model = fit2008_sex, dtms = definitiondtms,
                                    controls = list(time = definitiondtms$timescale, Gender = 0))
Tp2008_male     <- dtms_matrix(probs = probs2008_male, dtms = definitiondtms)
S2008_male      <- dtms_start(dtms = definitiondtms, data = Model2008_sex,
                              start_time = 50, variables = list(Gender = 0))
LE2008_male_sex <- dtms_expectancy(dtms = definitiondtms, matrix = Tp2008_male,
                                   start_distr = S2008_male, start_time = 50)
cat("\nLE50 2008-2018 Males (additive):\n"); print(LE2008_male_sex)

# 2008-2018 Females
probs2008_female  <- dtms_transitions(model = fit2008_sex, dtms = definitiondtms,
                                      controls = list(time = definitiondtms$timescale, Gender = 1))
Tp2008_female     <- dtms_matrix(probs = probs2008_female, dtms = definitiondtms)
S2008_female      <- dtms_start(dtms = definitiondtms, data = Model2008_sex,
                                start_time = 50, variables = list(Gender = 1))
LE2008_female_sex <- dtms_expectancy(dtms = definitiondtms, matrix = Tp2008_female,
                                     start_distr = S2008_female, start_time = 50)
cat("\nLE50 2008-2018 Females (additive):\n"); print(LE2008_female_sex)

# Bootstrap function for additive model (stacks Males/Females)
bootfun_sex <- function(data, dtms) {
  fit <- dtms_fit(data = data, controls = c("time", "Gender"))
  probs_m <- dtms_transitions(dtms = dtms, model = fit,
                              controls = list(time = dtms$timescale, Gender = 0))
  Tp_m    <- dtms_matrix(dtms = dtms, probs = probs_m)
  S_m     <- dtms_start(dtms = dtms, data = data, variables = list(Gender = 0))
  LE_m    <- dtms_expectancy(dtms = dtms, matrix = Tp_m, start_distr = S_m)
  probs_f <- dtms_transitions(dtms = dtms, model = fit,
                              controls = list(time = dtms$timescale, Gender = 1))
  Tp_f    <- dtms_matrix(dtms = dtms, probs = probs_f)
  S_f     <- dtms_start(dtms = dtms, data = data, variables = list(Gender = 1))
  LE_f    <- dtms_expectancy(dtms = dtms, matrix = Tp_f, start_distr = S_f)
  rbind(LE_m, LE_f)  # Males on top, Females on bottom
}

cat("\nBootstrapping 1996-2006 (additive, 1000 reps)...\n")
boot1996_sex <- dtms_boot(data = Model1996_sex, dtms = definitiondtms,
                          fun = bootfun_sex, idvar = "id", rep = 1000,
                          method = "block", parallel = TRUE)
print(summary(boot1996_sex))

cat("\nBootstrapping 2008-2018 (additive, 1000 reps)...\n")
boot2008_sex <- dtms_boot(data = Model2008_sex, dtms = definitiondtms,
                          fun = bootfun_sex, idvar = "id", rep = 1000,
                          method = "block", parallel = TRUE)
print(summary(boot2008_sex))


# --- Approach (3): Sex x Age interaction ---
# time_Gender = time * Gender (built manually because dtms_transitions does
# not parse formula syntax).
# For Males   (Gender = 0): time_Gender = 0 at all ages (no effect).
# For Females (Gender = 1): time_Gender = age (effect scales with age).

Model1996_sex$time_Gender <- Model1996_sex$time * Model1996_sex$Gender
Model2008_sex$time_Gender <- Model2008_sex$time * Model2008_sex$Gender

cat("\n=== DTMS: sex x age interaction ===\n")

fit1996_interact <- dtms_fit(data = Model1996_sex,
                             controls = c("time", "Gender", "time_Gender"))
fit2008_interact <- dtms_fit(data = Model2008_sex,
                             controls = c("time", "Gender", "time_Gender"))

# 1996-2006 Males
probs1996_male_int <- dtms_transitions(
  model = fit1996_interact, dtms = definitiondtms,
  controls = list(time = definitiondtms$timescale, Gender = 0,
                  time_Gender = definitiondtms$timescale * 0))
Tp1996_male_int <- dtms_matrix(probs = probs1996_male_int, dtms = definitiondtms)
S1996_male_int  <- dtms_start(dtms = definitiondtms, data = Model1996_sex,
                              start_time = 50, variables = list(Gender = 0))
LE1996_male_int <- dtms_expectancy(dtms = definitiondtms, matrix = Tp1996_male_int,
                                   start_distr = S1996_male_int, start_time = 50)
cat("\nLE50 1996-2006 Males (interaction):\n"); print(LE1996_male_int)

# 1996-2006 Females
probs1996_female_int <- dtms_transitions(
  model = fit1996_interact, dtms = definitiondtms,
  controls = list(time = definitiondtms$timescale, Gender = 1,
                  time_Gender = definitiondtms$timescale * 1))
Tp1996_female_int <- dtms_matrix(probs = probs1996_female_int, dtms = definitiondtms)
S1996_female_int  <- dtms_start(dtms = definitiondtms, data = Model1996_sex,
                                start_time = 50, variables = list(Gender = 1))
LE1996_female_int <- dtms_expectancy(dtms = definitiondtms, matrix = Tp1996_female_int,
                                     start_distr = S1996_female_int, start_time = 50)
cat("\nLE50 1996-2006 Females (interaction):\n"); print(LE1996_female_int)

# 2008-2018 Males
probs2008_male_int <- dtms_transitions(
  model = fit2008_interact, dtms = definitiondtms,
  controls = list(time = definitiondtms$timescale, Gender = 0,
                  time_Gender = definitiondtms$timescale * 0))
Tp2008_male_int <- dtms_matrix(probs = probs2008_male_int, dtms = definitiondtms)
S2008_male_int  <- dtms_start(dtms = definitiondtms, data = Model2008_sex,
                              start_time = 50, variables = list(Gender = 0))
LE2008_male_int <- dtms_expectancy(dtms = definitiondtms, matrix = Tp2008_male_int,
                                   start_distr = S2008_male_int, start_time = 50)
cat("\nLE50 2008-2018 Males (interaction):\n"); print(LE2008_male_int)

# 2008-2018 Females
probs2008_female_int <- dtms_transitions(
  model = fit2008_interact, dtms = definitiondtms,
  controls = list(time = definitiondtms$timescale, Gender = 1,
                  time_Gender = definitiondtms$timescale * 1))
Tp2008_female_int <- dtms_matrix(probs = probs2008_female_int, dtms = definitiondtms)
S2008_female_int  <- dtms_start(dtms = definitiondtms, data = Model2008_sex,
                                start_time = 50, variables = list(Gender = 1))
LE2008_female_int <- dtms_expectancy(dtms = definitiondtms, matrix = Tp2008_female_int,
                                     start_distr = S2008_female_int, start_time = 50)
cat("\nLE50 2008-2018 Females (interaction):\n"); print(LE2008_female_int)

# Bootstrap for interaction model
bootfun_interact <- function(data, dtms) {
  data$time_Gender <- data$time * data$Gender
  fit <- dtms_fit(data = data, controls = c("time", "Gender", "time_Gender"))
  probs_m <- dtms_transitions(dtms = dtms, model = fit,
                              controls = list(time = dtms$timescale, Gender = 0,
                                              time_Gender = dtms$timescale * 0))
  Tp_m <- dtms_matrix(dtms = dtms, probs = probs_m)
  S_m  <- dtms_start(dtms = dtms, data = data, variables = list(Gender = 0))
  LE_m <- dtms_expectancy(dtms = dtms, matrix = Tp_m, start_distr = S_m)
  probs_f <- dtms_transitions(dtms = dtms, model = fit,
                              controls = list(time = dtms$timescale, Gender = 1,
                                              time_Gender = dtms$timescale * 1))
  Tp_f <- dtms_matrix(dtms = dtms, probs = probs_f)
  S_f  <- dtms_start(dtms = dtms, data = data, variables = list(Gender = 1))
  LE_f <- dtms_expectancy(dtms = dtms, matrix = Tp_f, start_distr = S_f)
  rbind(LE_m, LE_f)
}

cat("\nBootstrapping 1996-2006 (interaction, 1000 reps)...\n")
boot1996_interact <- dtms_boot(data = Model1996_sex, dtms = definitiondtms,
                               fun = bootfun_interact, idvar = "id", rep = 1000,
                               method = "block", parallel = TRUE)
print(summary(boot1996_interact))

cat("\nBootstrapping 2008-2018 (interaction, 1000 reps)...\n")
boot2008_interact <- dtms_boot(data = Model2008_sex, dtms = definitiondtms,
                               fun = bootfun_interact, idvar = "id", rep = 1000,
                               method = "block", parallel = TRUE)
print(summary(boot2008_interact))


# --- Comparison across the three sex approaches ---
cat("\n=== STROKE HISTORY TOTAL LE50 ACROSS APPROACHES ===\n")

cat("\n--- 1996-2006 ---\n")
cat(sprintf(" Separate models:    Males %.2f   Females %.2f   Diff %.2f\n",
            LE1996Male["start:Stroke History_50", "TOTAL"],
            LE1996Female["start:Stroke History_50", "TOTAL"],
            LE1996Female["start:Stroke History_50", "TOTAL"] -
              LE1996Male["start:Stroke History_50", "TOTAL"]))
cat(sprintf(" Sex additive:       Males %.2f   Females %.2f   Diff %.2f\n",
            LE1996_male_sex["start:Stroke History_50", "TOTAL"],
            LE1996_female_sex["start:Stroke History_50", "TOTAL"],
            LE1996_female_sex["start:Stroke History_50", "TOTAL"] -
              LE1996_male_sex["start:Stroke History_50", "TOTAL"]))
cat(sprintf(" Sex x Age interact: Males %.2f   Females %.2f   Diff %.2f\n",
            LE1996_male_int["start:Stroke History_50", "TOTAL"],
            LE1996_female_int["start:Stroke History_50", "TOTAL"],
            LE1996_female_int["start:Stroke History_50", "TOTAL"] -
              LE1996_male_int["start:Stroke History_50", "TOTAL"]))

cat("\n--- 2008-2018 ---\n")
cat(sprintf(" Separate models:    Males %.2f   Females %.2f   Diff %.2f\n",
            LE2008Male["start:Stroke History_50", "TOTAL"],
            LE2008Female["start:Stroke History_50", "TOTAL"],
            LE2008Female["start:Stroke History_50", "TOTAL"] -
              LE2008Male["start:Stroke History_50", "TOTAL"]))
cat(sprintf(" Sex additive:       Males %.2f   Females %.2f   Diff %.2f\n",
            LE2008_male_sex["start:Stroke History_50", "TOTAL"],
            LE2008_female_sex["start:Stroke History_50", "TOTAL"],
            LE2008_female_sex["start:Stroke History_50", "TOTAL"] -
              LE2008_male_sex["start:Stroke History_50", "TOTAL"]))
cat(sprintf(" Sex x Age interact: Males %.2f   Females %.2f   Diff %.2f\n",
            LE2008_male_int["start:Stroke History_50", "TOTAL"],
            LE2008_female_int["start:Stroke History_50", "TOTAL"],
            LE2008_female_int["start:Stroke History_50", "TOTAL"] -
              LE2008_male_int["start:Stroke History_50", "TOTAL"]))

cat("\n--- GAINS (2008-2018 minus 1996-2006) ---\n")
cat(sprintf(" Separate models:    Males %.2f   Females %.2f\n",
            LE2008Male["start:Stroke History_50", "TOTAL"] -
              LE1996Male["start:Stroke History_50", "TOTAL"],
            LE2008Female["start:Stroke History_50", "TOTAL"] -
              LE1996Female["start:Stroke History_50", "TOTAL"]))
cat(sprintf(" Sex additive:       Males %.2f   Females %.2f\n",
            LE2008_male_sex["start:Stroke History_50", "TOTAL"] -
              LE1996_male_sex["start:Stroke History_50", "TOTAL"],
            LE2008_female_sex["start:Stroke History_50", "TOTAL"] -
              LE1996_female_sex["start:Stroke History_50", "TOTAL"]))
cat(sprintf(" Sex x Age interact: Males %.2f   Females %.2f\n",
            LE2008_male_int["start:Stroke History_50", "TOTAL"] -
              LE1996_male_int["start:Stroke History_50", "TOTAL"],
            LE2008_female_int["start:Stroke History_50", "TOTAL"] -
              LE1996_female_int["start:Stroke History_50", "TOTAL"]))


# ==============================================================================
# 18. DIAGNOSTICS: VIF, BMI, RACE/ETHNICITY (Tables S2, S3, S5)
# ==============================================================================

# --- Multicollinearity (Table S5) ---
# vif() doesn't accept interaction terms cleanly, so use simplified specs.
cat("\n=== MULTICOLLINEARITY DIAGNOSTICS (Table S5) ===\n")

cox_vif_male <- coxph(Surv(survival_timeYear, event) ~
                        stroke_baseline + stroke_acutecombo +
                        periodsince + Age + education_level +
                        smokev + hibpe + diabe +
                        heart + bmi_cat + insurance_type + race_eth,
                      data = filter(twosetofHRSedit, ragender == "Male"))

cox_vif_female <- coxph(Surv(survival_timeYear, event) ~
                          stroke_baseline + stroke_acutecombo +
                          periodsince + Age + education_level +
                          smokev + hibpe + diabe +
                          heart + bmi_cat + insurance_type + race_eth,
                        data = filter(twosetofHRSedit, ragender == "Female"))

cat("\nVIF - Males:\n");   print(vif(cox_vif_male))
cat("\nVIF - Females:\n"); print(vif(cox_vif_female))

cat("\nPairwise correlations between health variables:\n")
twosetofHRSedit %>%
  filter(!is.na(bmi_real)) %>%
  summarise(
    bmi_diabetes          = round(cor(bmi_real, diabe,             use = "complete.obs"), 3),
    bmi_hypertension      = round(cor(bmi_real, hibpe,             use = "complete.obs"), 3),
    bmi_heart             = round(cor(bmi_real, as.numeric(heart), use = "complete.obs"), 3),
    bmi_smoking           = round(cor(bmi_real, smokev,            use = "complete.obs"), 3),
    diabetes_hypertension = round(cor(diabe,    hibpe,             use = "complete.obs"), 3),
    diabetes_heart        = round(cor(diabe,    as.numeric(heart), use = "complete.obs"), 3),
    hypertension_heart    = round(cor(hibpe,    as.numeric(heart), use = "complete.obs"), 3)
  ) %>% print()

cat("\nBMI x Diabetes crosstab (row %):\n")
print(round(prop.table(table(twosetofHRSedit$bmi_cat, twosetofHRSedit$diabe),
                       margin = 1) * 100, 1))


# --- Mortality and demographics by race/ethnicity (Tables S2, S3) ---
cat("\n=== AGE-SPECIFIC MORTALITY BY RACE/ETHNICITY ===\n")
twosetofHRSedit %>%
  filter(!is.na(race_eth)) %>%
  mutate(age_group = cut(Age, breaks = c(50, 60, 70, 80, 90, 110))) %>%
  group_by(race_eth, age_group) %>%
  summarise(
    n                    = n(),
    n_deaths             = sum(event, na.rm = TRUE),
    mortality_per_1000py = round(n_deaths / sum(survival_timeYear, na.rm = TRUE) * 1000, 2),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = race_eth,
              values_from = c(n, n_deaths, mortality_per_1000py)) %>%
  print(width = Inf)

cat("\n=== AGE COMPOSITION BY RACE/ETHNICITY AND PERIOD ===\n")
twosetofHRSedit %>%
  filter(!is.na(race_eth)) %>%
  group_by(periodsince, race_eth) %>%
  summarise(
    n          = n_distinct(hhidpn),
    mean_age   = round(mean(Age, na.rm = TRUE), 1),
    sd_age     = round(sd(Age,   na.rm = TRUE), 1),
    pct_over75 = round(mean(Age > 75, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(periodsince, race_eth) %>%
  print()


# --- BMI mortality and age (Discussion: obesity paradox) ---
cat("\n=== MORTALITY RATES BY BMI CATEGORY ===\n")
twosetofHRSedit %>%
  filter(!is.na(bmi_cat)) %>%
  group_by(periodsince, ragender, bmi_cat) %>%
  summarise(
    n_obs                = n(),
    n_individuals        = n_distinct(hhidpn),
    n_deaths             = sum(event, na.rm = TRUE),
    total_py             = sum(survival_timeYear, na.rm = TRUE),
    mortality_per_1000py = round(n_deaths / total_py * 1000, 2),
    pct_died             = round(n_deaths / n_individuals * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(periodsince, ragender, bmi_cat) %>%
  print(n = 32)

cat("\n=== MEAN AGE BY BMI CATEGORY ===\n")
twosetofHRSedit %>%
  filter(!is.na(bmi_cat)) %>%
  group_by(periodsince, ragender, bmi_cat) %>%
  summarise(
    mean_age = round(mean(Age, na.rm = TRUE), 1),
    sd_age   = round(sd(Age,   na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(periodsince, ragender, bmi_cat) %>%
  print(n = 32)


# ==============================================================================
# 19. GAINS IN LE50 WITH 95% CIs (Supplementary Table S7)
# ==============================================================================
# Period-specific bootstraps are independent. SE of the gain combines them:
#    SE_gain = sqrt(SE_early^2 + SE_late^2),
#    SE_x    = (CI_hi - CI_lo) / (2 * 1.96)
# Hardcoded values are the point estimates and bootstrap CIs from the main
# analysis. Update them on re-runs.
# ==============================================================================

gain_ci <- function(le_early, ci_early_lo, ci_early_hi,
                    le_late,  ci_late_lo,  ci_late_hi) {
  se_early <- (ci_early_hi - ci_early_lo) / (2 * 1.96)
  se_late  <- (ci_late_hi  - ci_late_lo ) / (2 * 1.96)
  se_diff  <- sqrt(se_early^2 + se_late^2)
  gain     <- le_late - le_early
  lo       <- round(gain - 1.96 * se_diff, 2)
  hi       <- round(gain + 1.96 * se_diff, 2)
  cat(sprintf("  Gain: %.2f (%.2f to %.2f)\n", gain, lo, hi))
  c(gain = gain, lo = lo, hi = hi)
}

cat("\n=== GAINS IN LE50 WITH 95% CIs ===\n")

cat("\nBOTH SEXES\n")
cat(" No stroke start:\n");      g1 <- gain_ci(29.30, 29.01, 29.60, 30.42, 30.13, 30.71)
cat(" Stroke history start:\n"); g2 <- gain_ci(17.54, 16.61, 18.55, 20.59, 19.68, 21.67)
cat(" Average population:\n");   g3 <- gain_ci(29.02, 28.62, 29.33, 30.12, 29.81, 30.48)

cat("\nMALES\n")
cat(" No stroke start:\n");      g4 <- gain_ci(27.09, 26.77, 27.45, 28.77, 28.34, 29.23)
cat(" Stroke history start:\n"); g5 <- gain_ci(16.54, 15.43, 17.43, 19.92, 18.89, 20.87)
cat(" Average population:\n");   g6 <- gain_ci(26.79, 26.31, 27.19, 28.45, 28.06, 28.89)

cat("\nFEMALES\n")
cat(" No stroke start:\n");      g7 <- gain_ci(31.15, 30.85, 31.55, 31.76, 31.45, 32.09)
cat(" Stroke history start:\n"); g8 <- gain_ci(18.92, 17.82, 20.38, 21.44, 20.37, 22.72)
cat(" Average population:\n");   g9 <- gain_ci(30.90, 30.58, 31.33, 31.48, 31.12, 31.81)


# ==============================================================================
# 20. FIGURE 1: GAINS IN LE50 BY STARTING STATE AND SEX
# ==============================================================================

gains_data <- data.frame(
  group = rep(c("No Previous Stroke history",
                "Baseline Stroke History",
                "Average Population"), each = 3),
  sex   = rep(c("Both sexes", "Males", "Females"), 3),
  gain  = c(1.1, 1.7, 0.6,   3.1, 3.4, 2.5,   1.1, 1.7, 0.6),
  lo    = c(0.7, 1.1, 0.1,   1.7, 2.0, 0.8,   0.6, 1.1, 0.1),
  hi    = c(1.5, 2.2, 1.1,   4.4, 4.8, 4.3,   1.6, 2.3, 1.1)
)
gains_data$group <- factor(gains_data$group,
                           levels = c("No Previous Stroke history",
                                      "Baseline Stroke History",
                                      "Average Population"))
gains_data$sex   <- factor(gains_data$sex,
                           levels = c("Both sexes", "Males", "Females"))

sex_colors <- c("Both sexes" = "#6B6B6B", "Males" = "#2E5E8E", "Females" = "#C44E3E")

p <- ggplot(gains_data, aes(x = group, y = gain, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75),
           width = 0.65, color = "white", size = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                position = position_dodge(width = 0.75),
                width = 0.15, linewidth = 0.5, color = "#333333") +
  geom_text(aes(y = hi, label = sprintf("+%.1f", gain)),
            position = position_dodge(width = 0.75),
            vjust = -0.5, size = 4.5, fontface = "bold", color = "#333333") +
  scale_fill_manual(values = sex_colors) +
  scale_y_continuous(breaks = 0:5, labels = paste0("+", 0:5),
                     limits = c(0, 5.8), expand = c(0, 0)) +
  labs(y = "Gain in LE50 between 1996-2006 and 2008-2018 (years)",
       x = NULL, fill = NULL,
       caption = "Error bars = bootstrap 95% confidence intervals.") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    axis.line.x        = element_line(color = "grey40", linewidth = 0.3),
    axis.line.y        = element_line(color = "grey40", linewidth = 0.3),
    legend.position    = c(0.15, 0.92),
    legend.background  = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    legend.key.size    = unit(0.5, "cm"),
    legend.text        = element_text(size = 13),
    plot.caption       = element_text(size = 11, color = "grey50",
                                      hjust = 0.5, face = "italic"),
    axis.text          = element_text(size = 13),
    axis.title.y       = element_text(size = 14, face = "bold")
  )

# JAGS specs: PDF or EPS, 600 dpi, 2-column = 7.047 inches wide
pdf("Figure1_LE50_Gains.pdf", width = 7.047, height = 4.5)
print(p)
dev.off()

postscript("Figure1_LE50_Gains.eps", width = 7.047, height = 4.5,
           horizontal = FALSE, paper = "special")
print(p)
dev.off()


# ==============================================================================
# 21. EMPIRICAL CONFIRMATION OF HRS REFRESHMENT (Tables S10, S11)
# ==============================================================================
# (A) New entrants per HRS wave with mean entry age. Confirms steady-state
#     refreshment design (Sonnega et al., 2014).
# (B) Birth-cohort composition of age-50 starting observations per period.
#     Shows that 2008-2018 starting sample is largely independent of 1996.
# ==============================================================================

# --- (A) First-time entrants per wave (Table S10) ---
entry_info <- deathdatasetlong %>%
  filter(iwstat == 1) %>%
  group_by(hhidpn) %>%
  summarise(
    first_wave = min(wave),
    rabyear    = first(rabyear),
    .groups    = "drop"
  ) %>%
  mutate(age_at_entry = first_wave - rabyear)

cat("\n=== TABLE S10: First-time entrants by HRS wave ===\n")
entry_info %>%
  group_by(first_wave) %>%
  summarise(
    n_new_entrants = n(),
    mean_age_entry = round(mean(age_at_entry, na.rm = TRUE), 1),
    min_age_entry  = min(age_at_entry,  na.rm = TRUE),
    max_age_entry  = max(age_at_entry,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(first_wave) %>%
  print()

# --- (B) Birth-cohort composition of age-50 starting observations (Table S11) ---
deathdatasetlong <- deathdatasetlong %>%
  mutate(
    hrs_cohort = case_when(
      rabyear <= 1923                   ~ "Born 1923 or earlier",
      rabyear >= 1924 & rabyear <= 1930 ~ "Born 1924-1930",
      rabyear >= 1931 & rabyear <= 1941 ~ "Born 1931-1941",
      rabyear >= 1942 & rabyear <= 1947 ~ "Born 1942-1947",
      rabyear >= 1948 & rabyear <= 1953 ~ "Born 1948-1953",
      rabyear >= 1954 & rabyear <= 1959 ~ "Born 1954-1959",
      rabyear >= 1960 & rabyear <= 1965 ~ "Born 1960-1965",
      rabyear >= 1966                   ~ "Born 1966+",
      TRUE                              ~ NA_character_
    )
  )

age50_comp <- deathdatasetlong %>%
  filter(Age >= 50 & Age < 52) %>%
  filter(wave >= 1996 & wave <= 2018) %>%
  mutate(period = ifelse(wave <= 2006, "1996-2006", "2008-2018")) %>%
  distinct(hhidpn, period, .keep_all = TRUE)

cat("\n=== TABLE S11: Birth-cohort composition of age-50 starting observations ===\n")
tab <- table(age50_comp$hrs_cohort, age50_comp$period, useNA = "ifany")
print(tab)
cat("\nColumn percentages:\n")
print(round(prop.table(tab, margin = 2) * 100, 1))

# Headline statistic for response letter / Methods section
ids_in_1996_sample <- deathdatasetlong %>%
  filter(wave <= 2006 & iwstat == 1) %>%
  distinct(hhidpn) %>% pull(hhidpn)

age50_2008  <- age50_comp %>% filter(period == "2008-2018")
not_in_1996 <- sum(!age50_2008$hhidpn %in% ids_in_1996_sample)
cat(sprintf(
  "\nAge-50 obs in 2008-2018: %d, of which NOT in 1996 sample: %d (%.1f%%)\n",
  nrow(age50_2008), not_in_1996,
  100 * not_in_1996 / nrow(age50_2008)))


# ==============================================================================
# SESSION INFO (for OSF reproducibility)
# ==============================================================================
sessionInfo()
