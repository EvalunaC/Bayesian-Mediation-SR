suppressPackageStartupMessages({
  library(rjags)
  library(coda)
  library(survival)
  library(mvtnorm)
  library(EnvStats)
  library(HDInterval)
})

ARM_NAME <- "IMRT1Protons0"
MED_NAME <- "G4RIL"
SEED <- 2026
CAP_MONTHS <- 72
TSEQ <- seq(6, 72, 6)
LEVEL <- 0.95
MCMC_PRIMARY <- list(niter = 10000, nburn = 5000, nchain = 4, nthin = 5)#list(niter = 30000, nburn = 15000, nchain = 4, nthin = 15)
MCMC_SENS <- list(niter = 10000, nburn = 5000, nchain = 4, nthin = 5)
PEXP_KNOTS <- c(3, 6, 12, 24, 48)

gauss_legendre <- function(n) {
  i <- 1:(n - 1)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- b
  J[cbind(i + 1, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(nodes = e$values, weights = 2 * e$vectors[1, ]^2)
}
GL64 <- gauss_legendre(64)

pw_cuts <- function(knots) c(0, knots, Inf)

pw_exposure <- function(t, knots) {
  cuts <- pw_cuts(knots)
  lo <- cuts[-length(cuts)]
  hi <- cuts[-1]
  pmax(0, pmin(t, hi) - lo)
}

pw_interval <- function(t, knots) {
  cuts <- pw_cuts(knots)
  K <- length(cuts) - 1
  idx <- findInterval(t, cuts, rightmost.closed = FALSE, left.open = TRUE)
  idx[t <= 0] <- 1
  out <- matrix(0, length(t), K)
  out[cbind(seq_along(t), pmin(idx, K))] <- 1
  out
}

weib_loglik <- function(par, T1, d1, T2, d2, mat12, mat13, mat23) {
  p12 <- ncol(mat12); p13 <- ncol(mat13); p23 <- ncol(mat23)
  pp <- p12 + p13 + p23
  beta12 <- par[1:p12]
  beta13 <- par[(p12 + 1):(p12 + p13)]
  beta23 <- par[(p12 + p13 + 1):pp]
  nv12 <- exp(par[pp + 1]); nv13 <- exp(par[pp + 2]); nv23 <- exp(par[pp + 3])
  lamb12 <- exp(par[pp + 4]); lamb13 <- exp(par[pp + 5]); lamb23 <- exp(par[pp + 6])
  lp1 <- -lamb12 * T1^nv12 * exp(mat12 %*% beta12) - lamb13 * T1^nv13 * exp(mat13 %*% beta13)
  lp2 <- d2 * (1 - d1) * (log(nv13) + log(lamb13) + (nv13 - 1) * log(T2) + mat13 %*% beta13)
  lp3 <- d1 * (log(nv12) + log(lamb12) + (nv12 - 1) * log(T1) + mat12 %*% beta12 -
                 lamb23 * (T2^nv23 - T1^nv23) * exp(mat23 %*% beta23))
  lp4 <- d1 * d2 * (log(nv23) + log(lamb23) + (nv23 - 1) * log(T2) + mat23 %*% beta23)
  -sum(lp1 + lp2 + lp3 + lp4)
}

fit_mstate <- function(formu12, formu13, formu23, dat, model_file,
                       mcmc = MCMC_PRIMARY, knots = NULL, seed = SEED) {
  set.seed(seed)
  mat12 <- model.matrix(formu12, dat)[, -1, drop = FALSE]
  mat13 <- model.matrix(formu13, dat)[, -1, drop = FALSE]
  mat23 <- model.matrix(formu23, dat)[, -1, drop = FALSE]
  pfs <- model.extract(model.frame(formu12, dat), "response")
  os <- model.extract(model.frame(formu13, dat), "response")
  T1 <- as.numeric(pfs[, 1]); d1 <- as.numeric(pfs[, 2])
  T2 <- as.numeric(os[, 1]); d2 <- as.numeric(os[, 2])
  p12 <- ncol(mat12); p13 <- ncol(mat13); p23 <- ncol(mat23)
  pp <- p12 + p13 + p23
  N <- nrow(dat)
  nch <- mcmc$nchain

  fn_safe <- function(par, ...) {
    v <- weib_loglik(par, ...)
    if (is.finite(v)) v else 1e10
  }
  fit0 <- optim(rep(0, 6 + pp), fn_safe, control = list(maxit = 1e6),
                T1 = T1, d1 = d1, T2 = T2, d2 = d2,
                mat12 = mat12, mat13 = mat13, mat23 = mat23)
  par <- fit0$par
  hess <- try(optimHess(par, fn_safe, T1 = T1, d1 = d1, T2 = T2, d2 = d2,
                        mat12 = mat12, mat13 = mat13, mat23 = mat23), silent = TRUE)
  if (inherits(hess, "try-error") || any(!is.finite(hess))) hess <- diag(100, 6 + pp)
  vc <- function(idx) {
    v <- try(solve(hess[idx, idx, drop = FALSE]), silent = TRUE)
    if (inherits(v, "try-error")) diag(0.01, length(idx)) else v
  }
  b12i <- rmvnorm(nch, par[1:p12], vc(1:p12))
  b13i <- rmvnorm(nch, par[(p12 + 1):(p12 + p13)], vc((p12 + 1):(p12 + p13)))
  b23i <- rmvnorm(nch, par[(p12 + p13 + 1):pp], vc((p12 + p13 + 1):pp))

  datls <- list(N = N, p12 = p12, p13 = p13, p23 = p23,
                mat12 = mat12, mat13 = mat13, mat23 = mat23,
                T1 = T1, d1 = d1, T2 = T2, d2 = d2,
                C = 100, ones = rep(1, N))
  inits <- vector("list", nch)
  if (is.null(knots)) {
    for (ch in 1:nch) {
      inits[[ch]] <- list(beta12 = b12i[ch, ], beta13 = b13i[ch, ], beta23 = b23i[ch, ],
                          nv12 = exp(rnorm(1, par[pp + 1], 0.05)),
                          nv13 = exp(rnorm(1, par[pp + 2], 0.05)),
                          nv23 = exp(rnorm(1, par[pp + 3], 0.05)),
                          loglamb12 = rnorm(1, par[pp + 4], 0.1),
                          loglamb13 = rnorm(1, par[pp + 5], 0.1),
                          loglamb23 = rnorm(1, par[pp + 6], 0.1),
                          .RNG.name = "base::Mersenne-Twister", .RNG.seed = seed + ch)
    }
    pars <- c("nv12", "nv13", "nv23", "lamb12", "lamb13", "lamb23",
              "beta12", "beta13", "beta23")
    baseline <- list(type = "weibull", knots = NULL)
  } else {
    K <- length(knots) + 1
    datls$K <- K
    datls$E1 <- t(vapply(T1, pw_exposure, numeric(K), knots = knots))
    datls$W1 <- pw_interval(T1, knots)
    datls$W2 <- pw_interval(T2, knots)
    E2 <- t(vapply(T2, pw_exposure, numeric(K), knots = knots))
    datls$E23 <- pmax(E2 - datls$E1, 0)
    rate_init <- function(evt_time, evt_ind, expo) {
      ev <- colSums(pw_interval(evt_time, knots) * evt_ind)
      pmax(ev, 0.5) / pmax(colSums(expo), 1)
    }
    r12 <- rate_init(T1, d1, datls$E1)
    r13 <- rate_init(T2, d2 * (1 - d1), datls$E1)
    r23 <- rate_init(T2, d1 * d2, datls$E23)
    for (ch in 1:nch) {
      jit <- exp(rnorm(K, 0, 0.05))
      inits[[ch]] <- list(beta12 = rep(0, p12), beta13 = rep(0, p13), beta23 = rep(0, p23),
                          lam12v = r12 * jit, lam13v = r13 * jit, lam23v = r23 * jit,
                          .RNG.name = "base::Mersenne-Twister", .RNG.seed = seed + ch)
    }
    pars <- c("lam12v", "lam13v", "lam23v", "beta12", "beta13", "beta23")
    baseline <- list(type = "pexp", knots = knots)
  }

  m <- jags.model(model_file, data = datls, inits = inits,
                  n.chains = nch, n.adapt = 1000, quiet = TRUE)
  update(m, mcmc$nburn)
  ch <- coda.samples(m, pars, n.iter = mcmc$niter, thin = mcmc$nthin)

  S <- do.call(rbind, lapply(ch, as.matrix))
  summ_tr <- function(tr, p, mat) {
    b <- S[, paste0("beta", tr, "[", 1:p, "]"), drop = FALSE]
    hd <- t(apply(b, 2, function(x) hdi(x, credMass = LEVEL)))
    out <- cbind(Mean = colMeans(b), StdDev = apply(b, 2, sd),
                 HPD.lower = hd[, 1], HPD.upper = hd[, 2])
    rownames(out) <- paste0("Trans ", substr(tr, 1, 1), "->", substr(tr, 2, 2), ": ",
                            colnames(mat))
    out
  }
  list(data = dat, data.lst = datls,
       formus = list(formu12 = formu12, formu13 = formu13, formu23 = formu23),
       chains = ch,
       summary = list(summ12 = summ_tr("12", p12, mat12),
                      summ13 = summ_tr("13", p13, mat13),
                      summ23 = summ_tr("23", p23, mat23)),
       baseline = baseline, mcmc = mcmc, seed = seed)
}

rhs_formula <- function(formu) {
  reformulate(attr(delete.response(terms(formu)), "term.labels"))
}

get_coefs <- function(S, tr, p) {
  S[, paste0("beta", tr, "[", 1:p, "]"), drop = FALSE]
}

cov_weights <- function(fit, arm.name = ARM_NAME, med.name = MED_NAME) {
  dat <- fit$data
  mfs <- lapply(fit$formus, function(f) model.frame(rhs_formula(f), dat))
  covar.names <- unique(unlist(lapply(mfs, function(m) names(m))))
  covar.names <- covar.names[covar.names != med.name]
  w <- rep(1, nrow(dat))
  for (covnm in covar.names) {
    x <- dat[, covnm]
    numval <- length(levels(as.factor(x)))
    w <- w * demp(x, x, discrete = (numval < 10))
  }
  w
}

bsl_funs <- function(samples, baseline) {
  if (baseline$type == "weibull") {
    nv12 <- samples[, "nv12"]; nv13 <- samples[, "nv13"]; nv23 <- samples[, "nv23"]
    l12 <- samples[, "lamb12"]; l13 <- samples[, "lamb13"]; l23 <- samples[, "lamb23"]
    list(H12 = function(t) l12 * t^nv12, h12 = function(t) l12 * nv12 * t^(nv12 - 1),
         H13 = function(t) l13 * t^nv13, h13 = function(t) l13 * nv13 * t^(nv13 - 1),
         H23 = function(t) l23 * t^nv23)
  } else {
    K <- length(baseline$knots) + 1
    L12 <- samples[, paste0("lam12v[", 1:K, "]"), drop = FALSE]
    L13 <- samples[, paste0("lam13v[", 1:K, "]"), drop = FALSE]
    L23 <- samples[, paste0("lam23v[", 1:K, "]"), drop = FALSE]
    kn <- baseline$knots
    list(H12 = function(t) as.numeric(L12 %*% pw_exposure(t, kn)),
         h12 = function(t) L12[, which(pw_interval(t, kn)[1, ] == 1)],
         H13 = function(t) as.numeric(L13 %*% pw_exposure(t, kn)),
         h13 = function(t) L13[, which(pw_interval(t, kn)[1, ] == 1)],
         H23 = function(t) as.numeric(L23 %*% pw_exposure(t, kn)))
  }
}

lrr_mediation <- function(fit, tau, arm.name = ARM_NAME, med.name = MED_NAME,
                          gl = GL64, chunk = 1000) {
  dat <- fit$data
  f12 <- rhs_formula(fit$formus$formu12)
  f13 <- rhs_formula(fit$formus$formu13)
  f23 <- rhs_formula(fit$formus$formu23)
  mat12 <- model.matrix(f12, dat)[, -1, drop = FALSE]
  mat13 <- model.matrix(f13, dat)[, -1, drop = FALSE]
  mat23 <- model.matrix(f23, dat)[, -1, drop = FALSE]
  datc <- dat
  datc[, arm.name] <- 1
  mat12c <- model.matrix(f12, datc)[, -1, drop = FALSE]
  mat13c <- model.matrix(f13, datc)[, -1, drop = FALSE]
  mat23c <- model.matrix(f23, datc)[, -1, drop = FALSE]

  Arm <- dat[, arm.name]
  conind <- Arm == sort(unique(Arm))[1]
  trtind <- !conind
  w <- cov_weights(fit, arm.name, med.name)
  w1 <- w[trtind] / sum(w[trtind])
  w0 <- w[conind] / sum(w[conind])

  nchain <- length(fit$chains)
  perchain <- vector("list", nchain)
  for (cc in 1:nchain) {
    S_all <- as.matrix(fit$chains[[cc]])
    nS <- nrow(S_all)
    Med <- array(dim = c(length(tau), nS, 4))
    SP <- array(dim = c(length(tau), nS, 3))
    for (st in seq(1, nS, by = chunk)) {
      en <- min(st + chunk - 1, nS)
      rows <- st:en
      S <- S_all[rows, , drop = FALSE]
      p12 <- fit$data.lst$p12; p13 <- fit$data.lst$p13; p23 <- fit$data.lst$p23
      B12 <- get_coefs(S, "12", p12)
      B13 <- get_coefs(S, "13", p13)
      B23 <- get_coefs(S, "23", p23)
      bf <- bsl_funs(S, fit$baseline)
      ex12 <- exp(B12 %*% t(mat12)); ex13 <- exp(B13 %*% t(mat13)); ex23 <- exp(B23 %*% t(mat23))
      ex12c <- exp(B12 %*% t(mat12c)); ex13c <- exp(B13 %*% t(mat13c)); ex23c <- exp(B23 %*% t(mat23c))

      surv_at <- function(ti, e12, e13, e23) {
        u <- 0.5 * ti * (gl$nodes + 1)
        wq <- 0.5 * ti * gl$weights
        P11t <- exp(-bf$H12(ti) * e12 - bf$H13(ti) * e13)
        H23t <- bf$H23(ti)
        P01 <- 0
        for (g in seq_along(u)) {
          ug <- u[g]
          h12u <- bf$h12(ug)
          P00u <- exp(-bf$H12(ug) * e12 - bf$H13(ug) * e13)
          P22ut <- exp(-(H23t - bf$H23(ug)) * e23)
          P01 <- P01 + wq[g] * (P00u * (h12u * e12) * P22ut)
        }
        Sv <- P11t + P01
        Sv[Sv > 1] <- 1
        Sv
      }

      for (i in seq_along(tau)) {
        ti <- tau[i]
        Sf <- surv_at(ti, ex12, ex13, ex23)
        Sc <- surv_at(ti, ex12c, ex13c, ex23c)
        Survp1 <- as.numeric(Sf[, trtind, drop = FALSE] %*% w1)
        Survp0 <- as.numeric(Sf[, conind, drop = FALSE] %*% w0)
        Survp01 <- as.numeric(Sc[, conind, drop = FALSE] %*% w0)
        SP[i, rows, ] <- cbind(Survp0, Survp1, Survp01)
        Med[i, rows, ] <- cbind(log(Survp1) - log(Survp0),
                                log(Survp01) - log(Survp0),
                                log(Survp1) - log(Survp01),
                                (Survp1 - Survp01) / (Survp1 - Survp0))
      }
    }
    perchain[[cc]] <- list(Mediation = Med, SurvPred = SP)
  }
  bind2 <- function(...) {
    arrs <- list(...)
    d <- dim(arrs[[1]])
    out <- array(dim = c(d[1], sum(vapply(arrs, function(a) dim(a)[2], 0)), d[3]))
    at <- 0
    for (a in arrs) {
      out[, at + seq_len(dim(a)[2]), ] <- a
      at <- at + dim(a)[2]
    }
    out
  }
  Mediation <- do.call(bind2, lapply(perchain, `[[`, "Mediation"))
  SurvPred <- do.call(bind2, lapply(perchain, `[[`, "SurvPred"))
  list(Mediation = Mediation, SurvPred = SurvPred, perchain = perchain, tau = tau)
}

lrr_summary <- function(lrr, level = LEVEL) {
  lo <- (1 - level) / 2
  hi <- 1 - lo
  out <- do.call(rbind, lapply(seq_along(lrr$tau), function(i) {
    rows <- lapply(1:3, function(k) {
      x <- exp(lrr$Mediation[i, , k])
      data.frame(Time = lrr$tau[i], Effect = c("Total", "Direct", "Mediated")[k],
                 RR = median(x), RR_lo = quantile(x, lo), RR_hi = quantile(x, hi))
    })
    mp <- lrr$Mediation[i, , 4]
    rows[[4]] <- data.frame(Time = lrr$tau[i], Effect = "MedProp",
                            RR = median(mp, na.rm = TRUE),
                            RR_lo = quantile(mp, lo, na.rm = TRUE),
                            RR_hi = quantile(mp, hi, na.rm = TRUE))
    do.call(rbind, rows)
  }))
  rownames(out) <- NULL
  out
}

tranprob_med <- function(fit, tau, arm.name = ARM_NAME, med.name = MED_NAME,
                         gl = GL64, chunk = 1000, level = LEVEL) {
  dat <- fit$data
  f12 <- rhs_formula(fit$formus$formu12)
  f13 <- rhs_formula(fit$formus$formu13)
  f23 <- rhs_formula(fit$formus$formu23)
  mat12 <- model.matrix(f12, dat)[, -1, drop = FALSE]
  mat13 <- model.matrix(f13, dat)[, -1, drop = FALSE]
  mat23 <- model.matrix(f23, dat)[, -1, drop = FALSE]
  Arm <- dat[, arm.name]
  Medv <- dat[, med.name]
  grp <- list(IMRT_med0 = Arm == 1 & Medv == 0, IMRT_med1 = Arm == 1 & Medv == 1,
              PBT_med0 = Arm == 0 & Medv == 0, PBT_med1 = Arm == 0 & Medv == 1)
  w <- cov_weights(fit, arm.name, med.name)

  S_all <- do.call(rbind, lapply(fit$chains, as.matrix))
  nS <- nrow(S_all)
  qty <- c("RemainS", "StoP", "StoD", "PtoD")
  res <- array(dim = c(length(tau), nS, length(grp), length(qty)),
               dimnames = list(NULL, NULL, names(grp), qty))
  for (st in seq(1, nS, by = chunk)) {
    en <- min(st + chunk - 1, nS)
    rows <- st:en
    S <- S_all[rows, , drop = FALSE]
    p12 <- fit$data.lst$p12; p13 <- fit$data.lst$p13; p23 <- fit$data.lst$p23
    B12 <- get_coefs(S, "12", p12)
    B13 <- get_coefs(S, "13", p13)
    B23 <- get_coefs(S, "23", p23)
    bf <- bsl_funs(S, fit$baseline)
    ex12 <- exp(B12 %*% t(mat12)); ex13 <- exp(B13 %*% t(mat13)); ex23 <- exp(B23 %*% t(mat23))
    for (i in seq_along(tau)) {
      ti <- tau[i]
      u <- 0.5 * ti * (gl$nodes + 1)
      wq <- 0.5 * ti * gl$weights
      P11t <- exp(-bf$H12(ti) * ex12 - bf$H13(ti) * ex13)
      H23t <- bf$H23(ti)
      P12t <- 0
      P13t <- 0
      for (g in seq_along(u)) {
        ug <- u[g]
        P00u <- exp(-bf$H12(ug) * ex12 - bf$H13(ug) * ex13)
        P22ut <- exp(-(H23t - bf$H23(ug)) * ex23)
        P12t <- P12t + wq[g] * (P00u * (bf$h12(ug) * ex12) * P22ut)
        P13t <- P13t + wq[g] * (P00u * (bf$h13(ug) * ex13))
      }
      P23t <- 1 - exp(-H23t * ex23)
      for (gname in names(grp)) {
        gi <- grp[[gname]]
        wg <- w[gi] / sum(w[gi])
        res[i, rows, gname, "RemainS"] <- as.numeric(P11t[, gi, drop = FALSE] %*% wg)
        res[i, rows, gname, "StoP"] <- as.numeric(P12t[, gi, drop = FALSE] %*% wg)
        res[i, rows, gname, "StoD"] <- as.numeric(P13t[, gi, drop = FALSE] %*% wg)
        res[i, rows, gname, "PtoD"] <- as.numeric(P23t[, gi, drop = FALSE] %*% wg)
      }
    }
  }
  lo <- (1 - level) / 2
  hi <- 1 - lo
  summ <- expand.grid(Time = tau, Group = names(grp), Quantity = qty,
                      stringsAsFactors = FALSE)
  summ$Median <- summ$Lower <- summ$Upper <- NA_real_
  for (r in seq_len(nrow(summ))) {
    i <- match(summ$Time[r], tau)
    x <- res[i, , summ$Group[r], summ$Quantity[r]]
    summ$Median[r] <- median(x, na.rm = TRUE)
    summ$Lower[r] <- quantile(x, lo, na.rm = TRUE)
    summ$Upper[r] <- quantile(x, hi, na.rm = TRUE)
  }
  list(summary = summ, tau = tau, level = level)
}

diag_table <- function(fit) {
  ch <- fit$chains
  gd <- gelman.diag(ch, multivariate = FALSE, autoburnin = FALSE)
  ess <- effectiveSize(ch)
  data.frame(Parameter = rownames(gd$psrf),
             Rhat = gd$psrf[, 1], Rhat_upperCI = gd$psrf[, 2],
             ESS = as.numeric(ess[rownames(gd$psrf)]))
}

derived_diag <- function(lrr) {
  eff <- c("lRR_total", "lRR_direct", "lRR_mediated", "MedProp")
  clean <- function(x) x[is.finite(x)]
  out <- list()
  for (i in seq_along(lrr$tau)) {
    for (k in 1:4) {
      xs <- lapply(lrr$perchain, function(pc) clean(pc$Mediation[i, , k]))
      nmin <- min(vapply(xs, length, 0L))
      ml <- mcmc.list(lapply(xs, function(x) mcmc(x[seq_len(nmin)])))
      gd <- try(gelman.diag(ml, autoburnin = FALSE)$psrf[1, ], silent = TRUE)
      if (inherits(gd, "try-error")) gd <- c(NA, NA)
      out[[length(out) + 1]] <- data.frame(
        Time = lrr$tau[i], Quantity = sub("lRR_", "", eff[k]),
        Rhat = gd[1], Rhat_upperCI = gd[2],
        ESS = sum(vapply(xs, function(x) as.numeric(effectiveSize(mcmc(x))), 0)))
    }
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

evalue_rr <- function(rr) {
  r <- ifelse(rr < 1, 1 / rr, rr)
  r + sqrt(r * (r - 1))
}
