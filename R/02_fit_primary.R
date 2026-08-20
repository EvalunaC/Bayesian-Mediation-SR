source("R/00_functions.R")

dta <- readRDS("data/dta_matched_72.rds")
f12 <- Surv(Time_PFS, Recurrence_progression) ~ IMRT1Protons0 + G4RIL
f13 <- Surv(Time_OS, Status) ~ IMRT1Protons0 + G4RIL
f23 <- ~ IMRT1Protons0 + G4RIL

fit <- fit_mstate(f12, f13, f23, dta, "models/model_weibull.txt",
                  mcmc = MCMC_PRIMARY, seed = SEED)
saveRDS(fit, "output/fit_primary.rds")
