source("R/00_functions.R")

dta <- readRDS("data/dta_matched_72.rds")
zsc <- function(x) as.numeric(scale(x))
dta$zALC <- zsc(dta$CRT0ALC)
dta$zPTV <- zsc(dta$PTV)
dta$zChemo <- zsc(dta$Number_of_concurrent_chemotherapy_cycles)

mk <- function(rhs) list(
  as.formula(paste("Surv(Time_PFS, Recurrence_progression) ~", rhs)),
  as.formula(paste("Surv(Time_OS, Status) ~", rhs)),
  as.formula(paste("~", rhs)))

rt_start <- as.Date(dta$Radiation_Start_Date)
rt_end <- as.Date(dta$Radiation_End_Date)
rtdur <- as.numeric(rt_end - rt_start) / 30.4375
dlm <- dta
dlm$Time_PFS <- dta$Time_PFS - rtdur
dlm$Time_OS <- dta$Time_OS - rtdur
keep <- dlm$Time_PFS > 0 & dlm$Time_OS > 0
writeLines(c(
  sprintf("events_on_or_before_rt_end,%d",
          sum((dta$Recurrence_progression == 1 & dta$Time_PFS <= rtdur) |
              (dta$Status == 1 & dta$Time_OS <= rtdur))),
  sprintf("landmark_n_retained,%d", sum(keep))),
  "output/landmark_counts.csv")

specs <- list(
  prior_coef = list(dat = dta, rhs = "IMRT1Protons0 + G4RIL",
                    model = "models/model_weibull_prior_coef.txt", knots = NULL),
  prior_shape = list(dat = dta, rhs = "IMRT1Protons0 + G4RIL",
                     model = "models/model_weibull_prior_shape.txt", knots = NULL),
  pexp = list(dat = dta, rhs = "IMRT1Protons0 + G4RIL",
              model = "models/model_pexp.txt", knots = PEXP_KNOTS),
  ps = list(dat = dta, rhs = "IMRT1Protons0 + G4RIL + Propen.Score",
            model = "models/model_weibull.txt", knots = NULL),
  confounders = list(dat = dta, rhs = "IMRT1Protons0 + G4RIL + zALC + zPTV + zChemo",
                     model = "models/model_weibull.txt", knots = NULL),
  interaction = list(dat = dta, rhs = "IMRT1Protons0 * G4RIL",
                     model = "models/model_weibull.txt", knots = NULL),
  landmark = list(dat = dlm[keep, ], rhs = "IMRT1Protons0 + G4RIL",
                  model = "models/model_weibull.txt", knots = NULL))

tab <- list()
for (nm in names(specs)) {
  sp <- specs[[nm]]
  f <- mk(sp$rhs)
  fit <- fit_mstate(f[[1]], f[[2]], f[[3]], sp$dat, sp$model,
                    mcmc = MCMC_SENS, knots = sp$knots, seed = SEED)
  saveRDS(fit, sprintf("output/fit_sens_%s.rds", nm))
  dg <- diag_table(fit)
  lrr <- lrr_mediation(fit, TSEQ)
  saveRDS(lrr, sprintf("output/lrr_sens_%s.rds", nm))
  ss <- lrr_summary(lrr, level = LEVEL)
  ss$Spec <- nm
  ss$MaxRhat <- max(dg$Rhat, na.rm = TRUE)
  write.csv(ss, sprintf("output/mediation_sens_%s_95CrI.csv", nm), row.names = FALSE)
  tab[[nm]] <- subset(ss, Time == 72)
  if (nm == "interaction") {
    su <- rbind(fit$summary$summ12, fit$summary$summ13, fit$summary$summ23)
    write.csv(su[grepl("IMRT1Protons0:G4RIL", rownames(su)), , drop = FALSE],
              "output/interaction_coefficients.csv")
  }
}
write.csv(do.call(rbind, tab), "output/sensitivity_72mo_summary.csv", row.names = FALSE)
