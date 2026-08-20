# Code for "Bayesian Mediation Multi-State Model Elucidates the Impact of Cancer Immunity on Multiple Oncology Endpoints After Chemoradiotherapy for Esophageal Cancer" (Scientific Reports).

## Requirements

- R (>= 4.2)
- JAGS (>= 4.3)
- R packages: rjags, coda, survival, mvtnorm, EnvStats, HDInterval, ggplot2, gridExtra

## Data

The patient-level dataset contains protected health information and is not included. Place the propensity-matched cohort file at `data/df_propen_matched.csv` with columns: `Time_OS`, `Status`, `Time_PFS`, `Recurrence_progression`, `IMRT1Protons0` (1 = IMRT, 0 = PBT), `G4RIL`, `Propen.Score`, `CRT0ALC`, `PTV`, `Number_of_concurrent_chemotherapy_cycles`, `Radiation_Start_Date`, `Radiation_End_Date`.

## Run order

From this directory:

```
Rscript R/01_prepare_data.R
Rscript R/02_fit_primary.R
Rscript R/03_mediation_primary.R
Rscript R/04_diagnostics_gof.R
Rscript R/05_sensitivity.R
Rscript R/06_evalues.R
```

All outputs are written to `output/`.
