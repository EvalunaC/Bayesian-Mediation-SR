source("R/00_functions.R")

lrr <- readRDS("output/lrr_primary.rds")
ss <- lrr_summary(lrr, level = LEVEL)
ev <- subset(ss, Effect %in% c("Total", "Mediated") & Time %in% c(12, 36, 72))
ev$E_point <- evalue_rr(ev$RR)
ev$E_CI <- ifelse(ev$RR_hi >= 1 & ev$RR < 1, 1,
                  ifelse(ev$RR < 1, evalue_rr(ev$RR_hi), evalue_rr(ev$RR_lo)))
write.csv(ev, "output/evalues.csv", row.names = FALSE)
