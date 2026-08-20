source("R/00_functions.R")
suppressPackageStartupMessages({ library(ggplot2); library(gridExtra) })

fit <- readRDS("output/fit_primary.rds")

lrr <- lrr_mediation(fit, TSEQ)
saveRDS(lrr, "output/lrr_primary.rds")
ss <- lrr_summary(lrr, level = LEVEL)
write.csv(ss, "output/mediation_primary_95CrI.csv", row.names = FALSE)

pal <- c(Total = "#619CFF", Direct = "#F8766D", Mediated = "#00BA38")
eff <- subset(ss, Effect != "MedProp")
eff$Effect <- factor(eff$Effect, c("Total", "Direct", "Mediated"))
mp <- subset(ss, Effect == "MedProp")
p_rr <- ggplot(eff, aes(Time, RR, group = Effect, col = Effect, fill = Effect)) +
  geom_line(linewidth = 1.2) + geom_point() +
  geom_ribbon(aes(ymin = RR_lo, ymax = RR_hi), alpha = 0.2, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dotted", col = "black") +
  scale_color_manual(values = pal) + scale_fill_manual(values = pal) +
  xlab("Time in Months") + ylab("Risk Ratios") +
  theme_bw() + theme(legend.position = c(0.2, 0.25), legend.title = element_blank())
p_mp <- ggplot(mp, aes(Time, 100 * RR)) +
  geom_line(linewidth = 1.2, col = "grey40") + geom_point(col = "grey40") +
  scale_y_continuous(limits = c(0, 100)) +
  xlab("Time in Months") + ylab("Mediation Proportion (%)") + theme_bw()
ggsave("output/fig_mediation_primary.pdf",
       arrangeGrob(p_rr, p_mp, ncol = 2, widths = c(1.35, 1)), width = 10, height = 4.5)

tp <- tranprob_med(fit, tau = c(12, 24, 36, 48, 60, 72), level = LEVEL)
write.csv(tp$summary, "output/transprob_95CrI.csv", row.names = FALSE)
tpl <- tp$summary
tpl$Modality <- ifelse(grepl("IMRT", tpl$Group), "IMRT", "PBT")
tpl$G4RIL <- ifelse(grepl("med1", tpl$Group), "Yes", "No")
tpl$Quantity <- factor(tpl$Quantity, c("RemainS", "StoP", "StoD", "PtoD"),
                       c("Remain Stable", "Stable to Progression",
                         "Stable to Death", "Progression to Death"))
p_tp <- ggplot(tpl, aes(Time, Median, col = Modality, linetype = G4RIL)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.2) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Modality), alpha = 0.12, colour = NA) +
  facet_wrap(~Quantity, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c(IMRT = "#e96464", PBT = "#428ad2")) +
  scale_fill_manual(values = c(IMRT = "#e96464", PBT = "#428ad2")) +
  xlab("Time after start of CRT (months)") + ylab("Transition probability") +
  theme_bw() + theme(legend.position = "bottom")
ggsave("output/fig_transprob_primary.pdf", p_tp, width = 12, height = 3.6)

coefs <- rbind(
  data.frame(Transition = "S->P", fit$summary$summ12),
  data.frame(Transition = "S->D", fit$summary$summ13),
  data.frame(Transition = "P->D", fit$summary$summ23))
write.csv(coefs, "output/coefficients_primary.csv")
