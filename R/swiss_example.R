# swiss_example.R
#
# Data example for Section 6: the Swiss fertility and socioeconomic data of
# 1888 (R dataset `swiss`, from Mosteller and Tukey, 1977), 47 French-speaking
# provinces. The four variables given as percentages are treated as
# proportions on I = [0,1]. For each, the script reports the sample mean and
# variance (divisor n), the ratio S_n^2 / v_I(X̄_n), the statistic
# T_n = sqrt(n) G_n / omega(X̄_n), the asymptotic p-value 1 - Phi(T_n) of the
# test (14), a parametric-bootstrap p-value from 100,000 samples of size n
# from the truncated exponential with mean X̄_n, and the e-value (average
# likelihood ratio over 100 random splits) of the universal split LRT of
# Dunn et al., which rejects at level 0.05 when it is at least 20.
# Writes data/swiss_example.dat. Run from the repository root:
#   Rscript R/swiss_example.R

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(logcondens))

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
vI_omega2 <- function(m) {
  d <- Lder(Linv(1 - 2*m))
  c(d[2]/4, (d[4] + 2*d[2]^2 - d[3]^2/d[2])/16)
}
Tstat <- function(x) {
  n <- length(x); m <- mean(x); s2 <- mean((x - m)^2); vw <- vI_omega2(m)
  c(m, s2, s2/vw[1], sqrt(n) * (s2 - vw[1]) / sqrt(vw[2]))
}
# Vectorised T_n for a matrix of samples (one sample per row).
Tstat_rows <- function(X) {
  n <- ncol(X); m <- rowMeans(X); s2 <- rowMeans((X - m)^2)
  vw <- t(sapply(m, vI_omega2))
  sqrt(n) * (s2 - vw[, 1]) / sqrt(vw[, 2])
}
kde_reflect <- function(train, at) {
  h <- bw.nrd0(train)
  sapply(at, function(y) mean(dnorm((y - train)/h) + dnorm((y + train)/h) +
                              dnorm((y - (2 - train))/h)) / h)
}
e_value <- function(x, B = 100) {
  n <- length(x); r <- numeric(B)
  for (b in 1:B) {
    i0 <- sample.int(n, floor(n/2)); D0 <- x[i0]; D1 <- x[-i0]
    f0 <- logConDens(D0, smoothed = FALSE, print = FALSE)
    ll0 <- sum(evaluateLogConDens(D0, f0, which = 1)[, "log-density"])
    r[b] <- exp(min(sum(log(kde_reflect(D1, D0))) - ll0, 700))
  }
  mean(r)
}

con <- file(file.path("data", "swiss_example.dat"), "w")
writeLines("# variable | n | mean | variance | ratio | T_n | p asymptotic | p bootstrap | e-value", con)
for (nm in c("Agriculture", "Examination", "Education", "Catholic")) {
  x <- swiss[[nm]] / 100; n <- length(x)
  s <- Tstat(x)
  lam <- 2 * Linv(1 - 2*s[1])
  set.seed(1888)
  U <- matrix(runif(1e5 * n), ncol = n)
  Xb <- if (abs(lam) < 1e-12) U else -log1p(-U * (-expm1(-lam))) / lam
  pb <- mean(Tstat_rows(Xb) >= s[4])
  set.seed(1888)
  ev <- e_value(x)
  line <- sprintf("%-12s %d  %.3f  %.4f  %.3f  %6.2f  %.4f  %.4f  %.3g",
                  nm, n, s[1], s[2], s[3], s[4], 1 - pnorm(s[4]), pb, ev)
  writeLines(line, con); cat(line, "\n")
}
close(con)
