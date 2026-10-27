# Physical Activity, Adiposity & Hypertension

A reproducible **NHANES 2017–March 2020 pre-pandemic** analysis examining whether waist circumference helps explain the association between moderate-to-vigorous physical activity (MVPA) and hypertension.

![Analysis pathway](figures/analysis_pathway.svg)

## Research question

How is MVPA associated with hypertension, and how does the estimated MVPA association change after adding **waist circumference** to the outcome model?

The analysis uses a pathway-oriented, mediation-style regression sequence:

```
MVPA  ───────► Hypertension
  │                 ▲
  └──► Waist ───────┘
```

This repository presents the **analysis code and reproducible outputs only**. Classroom instructions, assignment rubrics, proposal documents, and student-facing guidance are intentionally excluded.

## Data

Public NHANES 2017–March 2020 pre-pandemic files are downloaded directly by the R script:

- Demographics
- Body measures
- Blood pressure
- Physical activity
- Blood-pressure questionnaire
- Smoking
- Alcohol use

The workflow focuses on adults ages **18–62**, uses complete cases for the modeled variables, and excludes participants reporting blood-pressure medication use.

## Variables

**Exposure:** weekly MVPA, combining moderate activity with vigorous activity weighted ×2.

**Mediator / pathway variable:** waist circumference (cm).

**Outcome:** hypertension indicator derived from repeated measured blood pressure, defined in the analysis as SBP ≥130 mmHg or DBP ≥80 mmHg.

**Covariates:** age, sex, race/ethnicity, family income-to-poverty ratio, sedentary time, smoking, and alcohol use.

## Modeling workflow

1. Download and harmonize seven NHANES components.
2. Construct weekly MVPA and sedentary time.
3. Average repeated systolic and diastolic BP measurements.
4. Build the complete-case analytic sample and exclude BP-medication users.
5. Fit pathway model H0 (linear): `waist ~ MVPA + covariates`.
6. Fit hypertension model H1: `HTN ~ MVPA + covariates`.
7. Fit hypertension model H2: `HTN ~ MVPA + waist + covariates`.
8. Compare H1 and H2 and describe how the MVPA coefficient changes when waist circumference enters the model; this is not an estimate of a mediated proportion.
9. Run regression diagnostics, VIF, ROC/AUC, Hosmer–Lemeshow checks, and influence diagnostics.

## Results

![Analytic sample results](figures/sample_results.svg)

The final analytic sample contains **247 adults** (**135 men and 112 women**). Mean age was **45.2 years**, mean waist circumference was **107.0 cm**, mean systolic blood pressure was **131.7 mmHg**, and **66.4%** met the analysis definition of hypertension. Mean weekly MVPA was **1,416.2 minute-equivalents**.

![Race and ethnicity distribution](figures/race_distribution.svg)

### Analysis snapshot

| Metric | Value |
|---|---:|
| Analytic sample | 247 |
| Mean age (SD) | 45.2 (11.3) years |
| Sex | 135 men; 112 women |
| Weekly MVPA (mean, SD) | 1,416 (1,974) min-equivalent |
| Waist circumference (mean, SD) | 107.0 (18.6) cm |
| Systolic BP (mean, SD) | 131.7 (20.8) mmHg |
| Hypertension prevalence | 66.4% |

These values describe the final complete-case, no-BP-medication analytic sample used in the project.

## Figures and generated outputs

The repository includes the descriptive result figures above. Running the analysis script also writes publication-ready model tables to `results/` and three portfolio-ready PNG figures to `figures/`:

- `mvpa_hypertension_models.png` — MVPA odds ratio before vs. after adding waist circumference.
- `mvpa_waist_association.png` — adjusted MVPA association with waist circumference.
- `sample_snapshot.png` — compact analytic-sample profile.

The coefficient figures are generated from the fitted models rather than hard-coded values.

## Repository structure

```
.
├── README.md
├── analysis/
│   └── nhanes_mvpa_hypertension.R
├── results/
│   └── sample_summary.csv
└── .gitignore
```

## Reproduce

Clone the repository, open R, install the required packages, and run:

```r
source("analysis/nhanes_mvpa_hypertension.R")
```

Required packages: `tidyverse`, `haven`, `DiagrammeR`, `broom`, `lmtest`, `car`, `pROC`, and `ResourceSelection`.

The script downloads the public NHANES XPT files directly from CDC/NCHS, so raw participant-level data are not stored in this repository.

## Limitations and interpretation

This is an **unweighted exploratory analysis** of a small selected complete-case sample, excluding participants taking blood-pressure medication. It is not nationally representative. Physical activity is self-reported; all variables are assessed cross-sectionally, and residual confounding or reverse causality remain possible. The three-model sequence evaluates whether the MVPA coefficient changes after waist circumference enters the hypertension model; logistic coefficient changes are also affected by non-collapsibility and should **not** be interpreted as a mediated share or formal causal mediation estimate. Numerical model coefficients should be read from generated outputs, not inferred from these descriptive summaries.

## Skills demonstrated

R · NHANES · multi-table data integration · feature engineering · epidemiologic modeling · linear regression · logistic regression · model diagnostics · reproducible analysis
