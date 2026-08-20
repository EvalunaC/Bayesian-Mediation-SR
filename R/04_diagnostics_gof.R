source("R/00_functions.R")
suppressPackageStartupMessages({ library(ggplot2) })

fit <- readRDS("output/fit_primary.rds")
lrr <- readRDS("output/lrr_primary.rds")
dta <- fit$data

write.csv(diag_table(fit), "output/diagnostics_parameters.csv", row.names = FALSE)
write.csv(derived_diag(lrr), "output/diagnostics_derived.csv", row.names = FALSE)

sel <- c("beta12[1]", "beta12[2]", "beta13[1]", "beta13[2]", "beta23[1]", "beta23[2]",
         "nv12", "nv13")
i72 <- match(72, lrr$tau)
mp_chains <- lapply(lrr$perchain, function(pc) pc$Mediation[i72, , 4])
pdf("output/fig_traces_primary.pdf", width = 14, height = 7.5)
op <- par(mfrow = c(3, 6), mar = c(2.4, 2.4, 2, 0.5), mgp = c(1.4, 0.4, 0))
cols <- c("#00000080", "#B2182B80", "#2166AC80", "#1B783780")
for (k in seq_along(sel)) {
  xs <- lapply(fit$chains, function(ch) as.numeric(ch[, sel[k]]))
  plot(xs[[1]], type = "n", ylim = range(unlist(xs)), xlab = "Iteration (thinned)",
       ylab = "", main = sel[k], cex.main = 0.8)
  for (cc in seq_along(xs)) lines(xs[[cc]], col = cols[cc], lwd = 0.5)
  dens <- lapply(xs, density)
  plot(NA, xlim = range(unlist(xs)), ylim = c(0, max(sapply(dens, function(d) max(d$y)))),
       xlab = "", ylab = "Density", main = "")
  for (cc in seq_along(dens)) lines(dens[[cc]], col = cols[cc], lwd = 1)
}
yl <- quantile(unlist(mp_chains), c(0.005, 0.995), na.rm = TRUE)
plot(mp_chains[[1]], type = "n", ylim = yl, xlab = "Iteration (thinned)", ylab = "",
     main = "Mediation proportion (72 mo)", cex.main = 0.8)
for (cc in seq_along(mp_chains)) lines(mp_chains[[cc]], col = cols[cc], lwd = 0.5)
dens <- lapply(mp_chains, function(x) density(x[is.finite(x)], from = yl[1], to = yl[2]))
plot(NA, xlim = yl, ylim = c(0, max(sapply(dens, function(d) max(d$y)))),
     xlab = "", ylab = "Density", main = "")
for (cc in seq_along(dens)) lines(dens[[cc]], col = cols[cc], lwd = 1)
par(op)
dev.off()

S <- do.call(rbind, lapply(fit$chains, as.matrix))
P <- list(nv12 = S[, "nv12"], nv13 = S[, "nv13"], nv23 = S[, "nv23"],
          l12 = S[, "lamb12"], l13 = S[, "lamb13"], l23 = S[, "lamb23"],
          b12a = S[, "beta12[1]"], b12m = S[, "beta12[2]"],
          b13a = S[, "beta13[1]"], b13m = S[, "beta13[2]"],
          b23a = S[, "beta23[1]"], b23m = S[, "beta23[2]"])
cell_curves <- function(a, m, tgrid, gl = GL64) {
  e12 <- exp(P$b12a * a + P$b12m * m)
  e13 <- exp(P$b13a * a + P$b13m * m)
  e23 <- exp(P$b23a * a + P$b23m * m)
  stab <- prog <- matrix(NA, length(P$nv12), length(tgrid))
  for (k in seq_along(tgrid)) {
    ti <- tgrid[k]
    if (ti == 0) { stab[, k] <- 1; prog[, k] <- 0; next }
    u <- 0.5 * ti * (gl$nodes + 1)
    wq <- 0.5 * ti * gl$weights
    H23t <- P$l23 * ti^P$nv23
    stab[, k] <- exp(-(P$l12 * ti^P$nv12 * e12 + P$l13 * ti^P$nv13 * e13))
    acc <- 0
    for (g in seq_along(u)) {
      ug <- u[g]
      h12u <- P$l12 * P$nv12 * ug^(P$nv12 - 1)
      P00u <- exp(-(P$l12 * ug^P$nv12 * e12 + P$l13 * ug^P$nv13 * e13))
      P22ut <- exp(-(H23t - P$l23 * ug^P$nv23) * e23)
      acc <- acc + wq[g] * (P00u * (h12u * e12) * P22ut)
    }
    prog[, k] <- acc
  }
  list(stable = stab, prog = prog, os = pmin(stab + prog, 1))
}

tg <- seq(0, 72, 2)
groups <- data.frame(a = c(1, 1, 0, 0), m = c(1, 0, 1, 0),
                     label = c("IMRT, G4RIL", "IMRT, non-G4RIL", "PBT, G4RIL", "PBT, non-G4RIL"))
curve_df <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
  cc <- cell_curves(groups$a[i], groups$m[i], tg)$os
  data.frame(t = tg, med = apply(cc, 2, median),
             lo = apply(cc, 2, quantile, 0.025), hi = apply(cc, 2, quantile, 0.975),
             Group = groups$label[i])
}))
km_df <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
  sub <- dta[dta$IMRT1Protons0 == groups$a[i] & dta$G4RIL == groups$m[i], ]
  sf <- survfit(Surv(Time_OS, Status) ~ 1, data = sub)
  data.frame(t = c(0, sf$time), s = c(1, sf$surv), Group = groups$label[i])
}))
p_gof <- ggplot() +
  geom_step(data = km_df, aes(t, s), colour = "black", linewidth = 0.4) +
  geom_ribbon(data = curve_df, aes(t, ymin = lo, ymax = hi), fill = "#4477AA", alpha = 0.20) +
  geom_line(data = curve_df, aes(t, med), colour = "#4477AA", linewidth = 0.9) +
  facet_wrap(~Group, nrow = 2) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Time after start of CRT (months)", y = "Overall survival probability") +
  theme_bw()
ggsave("output/fig_gof_os.pdf", p_gof, width = 9, height = 6.5)

mad_tab <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
  sub <- dta[dta$IMRT1Protons0 == groups$a[i] & dta$G4RIL == groups$m[i], ]
  sf <- survfit(Surv(Time_OS, Status) ~ 1, data = sub)
  kmf <- stepfun(sf$time, c(1, sf$surv))
  tt <- seq(2, min(72, max(sf$time)), 2)
  cc <- cell_curves(groups$a[i], groups$m[i], tt)$os
  data.frame(Group = groups$label[i], n = nrow(sub),
             MeanAbsDev = mean(abs(apply(cc, 2, median) - kmf(tt))))
}))
write.csv(mad_tab, "output/gof_mean_abs_dev.csv", row.names = FALSE)

wcell <- with(dta, table(IMRT1Protons0, G4RIL)) / nrow(dta)
occ <- Reduce(function(acc, i) {
  cc <- cell_curves(groups$a[i], groups$m[i], tg)
  w <- wcell[as.character(groups$a[i]), as.character(groups$m[i])]
  list(stable = acc$stable + w * cc$stable, prog = acc$prog + w * cc$prog)
}, seq_len(nrow(groups)), accumulate = FALSE, init = list(stable = 0, prog = 0))
occ_df <- rbind(
  data.frame(t = tg, med = apply(occ$stable, 2, median),
             lo = apply(occ$stable, 2, quantile, .025),
             hi = apply(occ$stable, 2, quantile, .975), State = "Stable"),
  data.frame(t = tg, med = apply(occ$prog, 2, median),
             lo = apply(occ$prog, 2, quantile, .025),
             hi = apply(occ$prog, 2, quantile, .975), State = "Progression (alive)"),
  data.frame(t = tg, med = apply(1 - occ$stable - occ$prog, 2, median),
             lo = apply(1 - occ$stable - occ$prog, 2, quantile, .025),
             hi = apply(1 - occ$stable - occ$prog, 2, quantile, .975), State = "Death"))
mk_long <- function(d) {
  T1 <- d$Time_PFS; d1 <- d$Recurrence_progression; T2 <- d$Time_OS; d2 <- d$Status
  T1a <- ifelse(d1 == 1 & T1 >= T2, pmax(T2 - 1e-6, 1e-6), T1)
  rows <- list()
  for (i in seq_len(nrow(d))) {
    if (d1[i] == 1) {
      rows[[length(rows) + 1]] <- data.frame(id = i, t0 = 0, t1 = T1a[i], ev = "prog")
      rows[[length(rows) + 1]] <- data.frame(id = i, t0 = T1a[i], t1 = T2[i],
                                             ev = ifelse(d2[i] == 1, "death", "censor"))
    } else {
      rows[[length(rows) + 1]] <- data.frame(id = i, t0 = 0, t1 = T2[i],
                                             ev = ifelse(d2[i] == 1, "death", "censor"))
    }
  }
  out <- do.call(rbind, rows)
  out$ev <- factor(out$ev, c("censor", "prog", "death"))
  out
}
long <- mk_long(dta)
sf <- survfit(Surv(t0, t1, ev) ~ 1, data = long, id = id)
aj <- data.frame(t = sf$time, sf$pstate)
names(aj) <- c("t", paste0("st_", sf$states))
aj_df <- rbind(
  data.frame(t = aj$t, p = aj[["st_(s0)"]], State = "Stable"),
  data.frame(t = aj$t, p = aj$st_prog, State = "Progression (alive)"),
  data.frame(t = aj$t, p = aj$st_death, State = "Death"))
slv <- c("Stable", "Progression (alive)", "Death")
occ_df$State <- factor(occ_df$State, slv)
aj_df$State <- factor(aj_df$State, slv)
p_occ <- ggplot() +
  geom_step(data = aj_df, aes(t, p), colour = "black", linewidth = 0.4) +
  geom_ribbon(data = occ_df, aes(t, ymin = lo, ymax = hi), fill = "#CC6677", alpha = 0.20) +
  geom_line(data = occ_df, aes(t, med), colour = "#CC6677", linewidth = 0.9) +
  facet_wrap(~State, nrow = 1) + coord_cartesian(xlim = c(0, 72), ylim = c(0, 1)) +
  labs(x = "Time after start of CRT (months)", y = "State-occupancy probability") +
  theme_bw()
ggsave("output/fig_occupancy.pdf", p_occ, width = 10, height = 3.6)

t13 <- ifelse(dta$Recurrence_progression == 1, dta$Time_PFS, dta$Time_OS)
e13 <- as.numeric(dta$Recurrence_progression == 0 & dta$Status == 1)
cz <- list(
  `1->2` = cox.zph(coxph(Surv(Time_PFS, Recurrence_progression) ~ IMRT1Protons0 + G4RIL, data = dta)),
  `1->3` = cox.zph(coxph(Surv(t13, e13) ~ IMRT1Protons0 + G4RIL, data = dta)),
  `2->3` = cox.zph(coxph(Surv(Time_OS - Time_PFS, Status) ~ IMRT1Protons0 + G4RIL,
                         data = dta[dta$Recurrence_progression == 1, ])))
ph_tab <- do.call(rbind, lapply(names(cz), function(nm) {
  tb <- cz[[nm]]$table
  data.frame(Transition = nm, Term = rownames(tb), ChiSq = tb[, "chisq"], p = tb[, "p"])
}))
write.csv(ph_tab, "output/ph_checks.csv", row.names = FALSE)
