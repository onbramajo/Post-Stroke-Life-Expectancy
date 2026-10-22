# Replication Code: Changes in Life Expectancy after Stroke among Americans Aged 50+, 1996–2018

This repository contains the R code to reproduce the main and supplementary analyses for:

> Bramajo O, Chakraborty M, Lisabeth L, Mehta NK. Changes in Life Expectancy after Stroke among Americans Aged 50+ by Sex, 1996–2018. *Journal of the American Geriatrics Society* (under review). Manuscript ID: JAGS-0093-CLI-Jan-26.

---

## Data

The analysis uses the **RAND HRS Longitudinal File (1992–2020)**. Data are not redistributed here due to the HRS data use agreement. Registered users can download the file at:

- HRS data portal: https://hrsdata.isr.umich.edu/ (search for "RAND HRS Longitudinal File")
- HRS registration: https://hrs.isr.umich.edu/data-products

After downloading, edit the path at the top of the script (`HRSdata <- read_dta("path/to/randhrs1992_2020v1.dta")`) to point to your local copy.

---

## Software

- R 4.x (tested on 4.3+)
- Required packages: `haven`, `tidyverse`, `dtms`, `survival`, `survminer`, `tableone`, `car`, `ggplot2`
- Full session info is printed at the end of the script via `sessionInfo()`

---

## Files

| File | Description |
|---|---|
| `HRS_code_main_analysis.R` | Main analysis: data preparation, Table 1, Cox proportional hazards models, DTMS life expectancy estimates (Sections 1–13). |
| `HRS_code_supplementary_sections_14_to_21.R` | Supplementary analyses: sensitivity tests, sex interactions, diagnostics, Figure 1, HRS refreshment confirmation (Sections 14–21). |

The two files are intended to run sequentially; the supplementary file uses objects created in the main file.

---

## Script structure

**Main analysis (Sections 1–13):**

1. Data preparation and variable construction
2. Reshape to long format
3. Age and survival time
4. Event indicator and recoding
5. Two cohorts: 1996–2006 and 2008–2018
6. Education categories
7. Stroke variables (baseline, incident, recurrent)
8. Data verification
9. Table 1 (descriptives)
10. Cox proportional hazards models (Tables 2 and 3)
11–12. Model comparison and fit statistics
13. DTMS life expectancy estimates (Table 4)

**Supplementary analyses (Sections 14–21):**

14. Random period assignment sensitivity (Table S8)
15. DTMS over full 1996–2018 period
16. Cox models with sex interactions (Table S9)
17. DTMS with sex as covariate / sex × age interaction (Table S11)
18. Multicollinearity, race/ethnicity, BMI diagnostics (Tables S2, S3, S5)
19. LE50 gains with 95% CIs (Table S7)
20. Figure 1
21. Empirical confirmation of HRS refreshment design (Tables S10, S11)

---

## Cox model numbering

The four Cox models in the paper map to these code objects:

| Paper | Code | Specification |
|---|---|---|
| Model 1 | `cox1_male`, `cox1_female` | Stroke history + period + age |
| Model 2 | `cox2_male`, `cox2_female` | + Incident stroke + interactions |
| Model 3 | `cox3_male`, `cox3_female` | + Sociodemographics, heart disease, BMI, insurance, race/ethnicity |
| Model 4 | `cox4_male`, `cox4_female` | + Period × covariate interactions |

---

## Run time

The full pipeline takes several hours on a standard laptop, mostly due to bootstrap iterations (`rep = 1000` in 13 separate `dtms_boot()` calls). To test that the script runs end-to-end without long waits, lower `rep` to 50 in the bootstrap calls.

---

## Reproducibility notes

- A fixed seed (`set.seed(12345)`) is used for the random period assignment sensitivity (Section 14).
- Bootstrap replications use a different random seed each run, so confidence interval bounds will vary slightly across runs. Point estimates are stable.

---

## Contact

Octavio Bramajo  
Institute for Evolutionary Medicine, University of Zurich  
octavio.bramajo@iem.uzh.ch

ORCID: [0000-0002-4148-6030](https://orcid.org/0000-0002-4148-6030)  
GitHub: [@onbramajo](https://github.com/onbramajo)  
OSF: [osf.io/e45hm](https://osf.io/e45hm/)
