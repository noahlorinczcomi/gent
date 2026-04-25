test_that("gent() returns a named list with the expected five components", {
  result <- gent(rnorm(5), diag(5))
  expect_named(result, c("pval", "shape", "rate", "mu_h0", "sigma2_h0"))
})

test_that("gent() returns length-1 numeric values for each component", {
  result <- gent(rnorm(5), diag(5))
  for (nm in names(result)) {
    expect_length(result[[nm]], 1L)
    expect_true(is.numeric(result[[nm]]))
  }
})

test_that("gent() p-value is in [0, 1]", {
  set.seed(1)
  result <- gent(rnorm(5), diag(5))
  expect_gte(result$pval, 0)
  expect_lte(result$pval, 1)
})

test_that("gent() null distribution parameters are correct for identity LD", {
  # With LD = I_m, tr(LD)=m, tr(LD^2)=m, so mu=m, sigma2=2m
  m <- 5
  result <- gent(rep(0, m), diag(m))
  expect_equal(result$mu_h0, m)
  expect_equal(result$sigma2_h0, 2 * m)
  # shape = mu^2 / sigma2 = m^2 / (2m) = m/2
  expect_equal(result$shape, m / 2)
  # rate = mu / sigma2 = m / (2m) = 1/2
  expect_equal(result$rate, 0.5)
})

test_that("gent() all-zero Z-statistics return p-value of 1", {
  result <- gent(rep(0, 5), diag(5))
  expect_equal(result$pval, 1)

  set.seed(7)
  R <- cov2cor(rWishart(1, 50, diag(5))[, , 1])
  result2 <- gent(rep(0, 5), R)
  expect_equal(result2$pval, 1)
})

test_that("gent() large Z-statistics yield a small p-value", {
  result <- gent(rep(5, 10), diag(10))
  expect_lt(result$pval, 0.05)
})

test_that("gent() shape and rate parameters are positive", {
  set.seed(42)
  result <- gent(rnorm(8), diag(8))
  expect_gt(result$shape, 0)
  expect_gt(result$rate, 0)
})

test_that("gent() with weight matrix A = I is identical to standard call", {
  set.seed(3)
  zs <- rnorm(5)
  r1 <- gent(zs, diag(5))
  r2 <- gent(zs, diag(5), A = diag(5))
  expect_equal(r1$pval,    r2$pval,    tolerance = 1e-10)
  expect_equal(r1$mu_h0,   r2$mu_h0,   tolerance = 1e-10)
  expect_equal(r1$sigma2_h0, r2$sigma2_h0, tolerance = 1e-10)
})

test_that("gent() with chisquares gives same p-value as passing zs directly", {
  zs <- c(2, 1, 0, 0, 0)
  r1 <- gent(zs, diag(5))
  r2 <- gent(zs, diag(5), chisquares = zs^2)
  expect_equal(r1$pval, r2$pval)
})

test_that("gent() returns a finite p-value for a correlated LD matrix", {
  set.seed(7)
  R <- cov2cor(rWishart(1, 50, diag(5))[, , 1])
  result <- gent(rnorm(5), R)
  expect_true(is.finite(result$pval))
})

test_that("gent() p-value decreases monotonically as signal strength grows", {
  LD <- diag(5)
  pvals <- sapply(c(0, 1, 3, 5, 8), function(mu) gent(rep(mu, 5), LD)$pval)
  expect_true(all(diff(pvals) <= 0))
})
