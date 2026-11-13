# power_comparison.R
#
# Compares, on I = [0,1], the diagnostic test (14) of Section 6 with the
# universal-inference split likelihood-ratio test for log-concavity of Dunn,
# Gangrade, Wasserman and Ramdas (JCGS, 2025), which is valid in finite
# samples. For the latter, the sample is split at random into halves D0 and
# D1; the numerator is a kernel density estimate from D1 (Gaussian kernel,
# reflected at 0 and 1, bandwidth bw.nrd0), the denominator is the
# log-concave maximum likelihood estimate from D0 (logcondens), both
# evaluated on D0, and the test rejects when the average of the likelihood
# ratio over B = 20 random splits is at least 1/alpha.
#
# alpha = 0.05, 2,000 replications per cell. Writes
# data/power_comparison.dat. Requires the R package logcondens.
# Run from the repository root:  Rscript R/power_comparison.R

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(logcondens))
library(parallel)

alpha <- 0.05
NS <- c(20, 50, 100, 500)
REPS <- 2000
B <- 20

# ---------------------------------------------------------------- Langevin
# L and its first three derivatives, with Taylor series near 0.
Lder <- function(z) {
  if (abs(z) < 0.05) {
    return(c(z/3 - z^3/45 + 2*z^5/945,
             1/3 - z^2/15 + 2*z^4/189,
             -2*z/15 + 8*z^3/189,
             -2/15 + 24*z^2/189))
  }
  c <- 1/tanh(z); q <- 1/sinh(z)^2
  c(c - 1/z, 1/z^2 - q, -2/z^3 + 2*q*c, 6/z^4 - 4*q*c^2 - 2*q^2)
}
Linv <- function(y) uniroot(function(z) Lder(z)[1] - y, c(-1e4, 1e4), tol = 1e-13)$root
# v_I(m) and omega(m)^2 = v_I^2 (2 + v_I'') on [0,1] (Remark 4).
vI_omega2 <- function(m) {
  d <- Lder(Linv(1 - 2*m))
  c(d[2]/4, (d[4] + 2*d[2]^2 - d[3]^2/d[2])/16)
}

boundary_test <- function(x) {
  n <- length(x); m <- mean(x); s2 <- mean((x - m)^2)
  vw <- vI_omega2(m)
  sqrt(n) * (s2 - vw[1]) > qnorm(1 - alpha) * sqrt(vw[2])
}

# ------------------------------------------------- universal split LRT
kde_reflect <- function(train, at) {
  h <- bw.nrd0(train)
  sapply(at, function(y) mean(dnorm((y - train)/h) + dnorm((y + train)/h) +
                              dnorm((y - (2 - train))/h)) / h)
}
universal_test <- function(x) {
  n <- length(x); ratio <- numeric(B)
  for (b in 1:B) {
    i0 <- sample.int(n, floor(n/2)); D0 <- x[i0]; D1 <- x[-i0]
    f0 <- logConDens(D0, smoothed = FALSE, print = FALSE)
    ll0 <- sum(evaluateLogConDens(D0, f0, which = 1)[, "log-density"])
    ll1 <- sum(log(kde_reflect(D1, D0)))
    ratio[b] <- exp(min(ll1 - ll0, 700))
  }
  mean(ratio) >= 1/alpha
}

# ---------------------------------------------------------------- densities
texp_sampler <- function(mu) {
  lam <- 2 * Linv(1 - 2*mu)
  if (abs(lam) < 1e-12) return(function(n) runif(n))
  function(n) -log1p(-runif(n) * (-expm1(-lam))) / lam
}
tnorm_sampler <- function(mu, sd) {
  a <- pnorm((0 - mu)/sd); b <- pnorm((1 - mu)/sd)
  function(n) mu + sd * qnorm(a + runif(n) * (b - a))
}
mix_beta <- function(a1, b1, a2, b2, w = 0.5) function(n) {
  k <- runif(n) < w
  ifelse(k, rbeta(n, a1, b1), rbeta(n, a2, b2))
}
dens <- list(
  list("truncated exponential, mean 0.25", texp_sampler(0.25)),
  list("uniform", function(n) runif(n)),
  list("truncated normal N(-2,1)", tnorm_sampler(-2, 1)),
  list("Beta(0.5,0.5)", function(n) rbeta(n, 0.5, 0.5)),
  list("Beta(0.8,0.8)", function(n) rbeta(n, 0.8, 0.8)),
  list("Beta(0.5,1.5)", function(n) rbeta(n, 0.5, 1.5)),
  list("mixture B(3,12),B(12,3)", mix_beta(3, 12, 12, 3)),
  list("mixture B(5,15),B(15,5)", mix_beta(5, 15, 15, 5)),
  list("mixture 0.8 B(2,8), 0.2 B(30,3)", mix_beta(2, 8, 30, 3, 0.8)))

ratio_of <- function(sampler) {
  set.seed(1); x <- sampler(2e6); m <- mean(x)
  var(x) / vI_omega2(m)[1]
}

run_cell <- function(k) {
  d <- dens[[k]]; out <- c()
  for (n in NS) {
    set.seed(20261113 + 100*k + n)
    rb <- ru <- 0
    for (r in 1:REPS) {
      x <- d[[2]](n)
      rb <- rb + boundary_test(x)
      ru <- ru + universal_test(x)
    }
    out <- c(out, rb/REPS, ru/REPS)
  }
  out
}

RNGkind("L'Ecuyer-CMRG")
res <- mclapply(seq_along(dens), run_cell, mc.cores = max(1, detectCores() - 1))
path <- file.path("data", "power_comparison.dat")
con <- file(path, "w")
writeLines(sprintf("# I = [0,1], alpha = %.2f, %d replications: density | ratio | (boundary, universal) for n = %s",
                   alpha, REPS, paste(NS, collapse = ", ")), con)
for (k in seq_along(dens)) {
  line <- sprintf("%-34s %.3f  %s", dens[[k]][[1]], ratio_of(dens[[k]][[2]]),
                  paste(sprintf("%.3f", res[[k]]), collapse = "  "))
  writeLines(line, con); cat(line, "\n")
}
close(con)
