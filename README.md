# CTMC-CKD

This repository contains the code, processed data, and results for a six-state continuous-time Markov chain (CTMC) model of Chronic Kidney Disease (CKD) progression.

The model uses CKD stages 1–5 plus death, with transitions only to the next stage or directly to death. Transition rates are fitted using longitudinal synthetic health-record data and a proportional hazards style parameterization.

## Project overview

- Disease model: 6-state CTMC with states 1–5 as CKD stages and state 6 as death.
- Assumptions:
  - CKD progression is irreversible.
  - Stage-skipping is not allowed.
  - Death is an absorbing state.
- Data included:
  - Processed cohort and panel files in `data_processed/`.
  - Result tables and plots in `results/`.
- Code included:
  - Data preparation and CTMC model scripts in `scripts/`.

## Repository structure

- `data_processed/`
  - Processed input data used by the analysis scripts.
  - Contains `ckd_panel.csv` and `cohort_summary.csv`.
- `results/`
  - Generated output from model fitting and evaluation.
  - Includes parameter tables, validation metrics, and plots.
- `scripts/`
  - Scripts for preprocessing, model fitting, cross-validation, and sensitivity analysis:
    - `scripts/preprocess.py`
    - `scripts/ckd_msm_M1.R`
    - `scripts/ckd_msm_M1_cv.R`
    - `scripts/ckd_msm_M1_sensitivity.R`
    - `scripts/run_ckd_msm_M1_cv.sh`
    - `scripts/run_ckd_msm_M1_cv_small.sh`

## Data generation

- Synthetic data are generated using Synthea.
- If you want to regenerate the raw data yourself, use the Synthea project at https://github.com/synthetichealth/synthea.
- After generating raw Synthea output, run `python3 scripts/preprocess.py` to produce the processed files used by this repository.
- If you do not want to generate your own data, use the already processed data in `data_processed/`.

## How to use

1. Clone the repository:
```bash
git clone https://github.com/ShahriyariLab/Multi-state-Continuous-Time-Markov-Chain-Modeling-for-Chronic-Kidney-Disease-Progression.git
cd Multi-state-Continuous-Time-Markov-Chain-Modeling-for-Chronic-Kidney-Disease-Progression
```
2. Inspect the processed data in `data_processed/`.
3. Run preprocessing if you need to regenerate inputs:
```bash
python3 scripts/preprocess.py
```
4. Fit the main CTMC model with R:
```bash
Rscript scripts/ckd_msm_M1.R
```
5. Run cross-validation:
```bash
Rscript scripts/ckd_msm_M1_cv.R
```
6. Run sensitivity analysis:
```bash
Rscript scripts/ckd_msm_M1_sensitivity.R
```

## Notes

- This repository only tracks `data_processed/`, `results/`, and `scripts/`.
- The manuscript source file `BMB/BMB.tex` is not uploaded to this repository.
- The `results/` folder may contain additional output files produced by running the scripts.
