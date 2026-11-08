# Pólya-Gamma sampler (src/include/stochtree/polya_gamma_sampler.h).
#
# Closed forms (Polson, Scott and Windle 2013):
#   E[PG(b, z)]   = b tanh(z / 2) / (2 z)
#   Var[PG(b, z)] = b (sinh(z) - z) / (4 z^3 cosh^2(z / 2))
#   E[exp(-t PG(b, z))] = cosh^b(z / 2) / cosh^b(sqrt((z^2 / 2 + t) / 2))
# The cumulants follow from the gamma-series representation:
#   kappa_n = b (n - 1)! (2 pi^2)^-n sum_k ((k - 1/2)^2 + z^2 / (4 pi^2))^-n.
#
# Tolerances. With N draws the sample mean has standard error sqrt(kappa_2 / N)
# and the sample variance sqrt((kappa_4 + 2 kappa_2^2) / N) to first order. We
# accept within 5 standard errors: each check then fails by chance with
# probability about 6e-7, so the whole file (about 50 checks) has a false-alarm
# rate near 3e-5 while still detecting a 2% bias in the variance at N = 20,000
# for every case below. The seed is fixed, so the tests are deterministic.

pg_mean <- function(b, z) {
  ifelse(abs(z) < 1e-6, b / 4, b * tanh(z / 2) / (2 * z))
}
pg_var <- function(b, z) {
  ifelse(
    abs(z) < 1e-6,
    b / 24,
    b * (sinh(z) - z) / (4 * z^3 * cosh(z / 2)^2)
  )
}
pg_cumulant <- function(b, z, n) {
  k <- seq_len(2e5)
  d <- (k - 0.5)^2 + z^2 / (4 * pi^2)
  b * factorial(n - 1) * (2 * pi^2)^(-n) * sum(d^(-n))
}
pg_laplace <- function(b, z, t) {
  (cosh(z / 2) / cosh(sqrt((z^2 / 2 + t) / 2)))^b
}

pg_cases <- list(
  c(b = 1, z = 0),
  c(b = 1, z = 2),
  c(b = 1, z = -8),
  c(b = 1, z = 60),
  c(b = 3, z = 0.5),
  c(b = 0.3, z = 0),
  c(b = 0.3, z = 4),
  c(b = 2.5, z = 1),
  c(b = 2.5, z = -20),
  c(b = 50, z = 0),
  c(b = 50, z = 7),
  c(b = 73.4, z = -2),
  c(b = 500, z = 0),
  c(b = 500, z = 3),
  c(b = 500, z = -40),
  c(b = 20000, z = 0),
  c(b = 20000, z = 5)
)

test_that("Polya-Gamma draws match the closed-form mean, variance and Laplace transform", {
  skip_on_cran()
  rng <- createCppRNG(20261108)
  n_draws <- 20000
  for (case in pg_cases) {
    b <- case[["b"]]
    z <- case[["z"]]
    x <- nashtree:::sample_polya_gamma_cpp(rep(b, n_draws), z, rng$rng_ptr)
    label <- sprintf("PG(%g, %g)", b, z)
    expect_true(all(is.finite(x) & x > 0), label = label)

    m <- pg_mean(b, z)
    v <- pg_var(b, z)
    k4 <- pg_cumulant(b, z, 4)
    # The closed form and the series agree; guards the reference itself
    expect_equal(pg_cumulant(b, z, 2), v, tolerance = 1e-4)

    se_mean <- sqrt(v / n_draws)
    se_var <- sqrt((k4 + 2 * v^2) / n_draws)
    expect_lt(abs(mean(x) - m), 5 * se_mean, label = paste(label, "mean"))
    expect_lt(abs(var(x) - v), 5 * se_var, label = paste(label, "variance"))

    # Laplace transform at t chosen so that t * E[x] is about 1
    t <- 1 / m
    lt <- exp(-t * x)
    expect_lt(
      abs(mean(lt) - pg_laplace(b, z, t)),
      5 * sd(lt) / sqrt(n_draws),
      label = paste(label, "Laplace transform")
    )
  }
})

test_that("Exact (integer b) and gamma-series (non-integer b) paths agree in distribution", {
  skip_on_cran()
  rng <- createCppRNG(7)
  n_draws <- 20000
  for (z in c(0, 1.5, -6)) {
    # b = 2 uses the exact sum of two Devroye draws; b = 2 + 1e-9 uses the series
    exact <- nashtree:::sample_polya_gamma_cpp(rep(2, n_draws), z, rng$rng_ptr)
    series <- nashtree:::sample_polya_gamma_cpp(rep(2 + 1e-9, n_draws), z, rng$rng_ptr)
    expect_gt(suppressWarnings(ks.test(exact, series)$p.value), 1e-4)
  }
})

test_that("sample_polya_gamma_cpp recycles, handles b = 0 and validates input", {
  skip_on_cran()
  rng <- createCppRNG(3)
  pg <- nashtree:::sample_polya_gamma_cpp
  expect_length(pg(c(1, 2, 3), 0.5, rng$rng_ptr), 3)
  expect_length(pg(2, c(-1, 0, 1, 2), rng$rng_ptr), 4)
  expect_length(pg(numeric(0), 1, rng$rng_ptr), 0)
  expect_equal(pg(0, 1, rng$rng_ptr), 0)
  expect_error(pg(c(1, 2), c(1, 2, 3), rng$rng_ptr), "same length")
  expect_error(pg(-1, 0, rng$rng_ptr), "non-negative")
  expect_error(pg(1, NA_real_, rng$rng_ptr), "finite")
  expect_error(pg(Inf, 0, rng$rng_ptr), "finite")

  # Same seed, same draws
  a <- pg(rep(c(1, 2.7, 600), 10), 1.3, createCppRNG(11)$rng_ptr)
  b <- pg(rep(c(1, 2.7, 600), 10), 1.3, createCppRNG(11)$rng_ptr)
  expect_identical(a, b)
})
