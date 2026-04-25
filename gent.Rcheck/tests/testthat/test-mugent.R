make_ldlist <- function(p, m, seed = 1) {
  set.seed(seed)
  lapply(seq_len(p), function(i) cov2cor(rWishart(1, m + 10, diag(m))[, , 1]))
}

make_Z <- function(m, p, seed = 1) {
  set.seed(seed)
  matrix(rnorm(m * p), nrow = m, ncol = p)
}

# ── mugent ────────────────────────────────────────────────────────────────────

test_that("mugent() returns a named list with the expected five components", {
  result <- mugent(make_Z(5, 2), make_ldlist(2, 5))
  expect_named(result, c("pval", "shape", "rate", "mu_h0", "sigma2_h0"))
})

test_that("mugent() p-value is in [0, 1]", {
  result <- mugent(make_Z(5, 2), make_ldlist(2, 5))
  expect_gte(result$pval, 0)
  expect_lte(result$pval, 1)
})

test_that("mugent() returns length-1 numeric components", {
  result <- mugent(make_Z(5, 2), make_ldlist(2, 5))
  for (nm in names(result)) {
    expect_length(result[[nm]], 1L)
    expect_true(is.numeric(result[[nm]]))
  }
})

test_that("mugent() all-zero Z-statistics return p-value of 1", {
  Z <- matrix(0, nrow = 5, ncol = 2)
  result <- mugent(Z, make_ldlist(2, 5))
  expect_equal(result$pval, 1)
})

test_that("mugent() large Z-statistics yield a small p-value", {
  Z <- matrix(5, nrow = 10, ncol = 2)
  result <- mugent(Z, make_ldlist(2, 10))
  expect_lt(result$pval, 0.05)
})

test_that("mugent() accepts a single (non-list) LD matrix", {
  set.seed(1)
  LD <- cov2cor(rWishart(1, 20, diag(5))[, , 1])
  result <- mugent(make_Z(5, 2), LD)
  expect_named(result, c("pval", "shape", "rate", "mu_h0", "sigma2_h0"))
})

test_that("mugent() mu_h0 equals number of populations p", {
  # Under H0, EZ = t(j) %*% I_p %*% j = sum(j) = p
  p <- 3
  result <- mugent(matrix(0, 5, p), make_ldlist(p, 5))
  expect_equal(result$mu_h0, p)
})

# ── mugent_ph ─────────────────────────────────────────────────────────────────

test_that("mugent_ph() returns a named list with the expected five components", {
  result <- mugent_ph(make_Z(5, 2), make_ldlist(2, 5))
  expect_named(result, c("pval", "shape", "rate", "mu_h0", "sigma2_h0"))
})

test_that("mugent_ph() p-value is in [0, 1]", {
  result <- mugent_ph(make_Z(5, 2), make_ldlist(2, 5))
  expect_gte(result$pval, 0)
  expect_lte(result$pval, 1)
})

test_that("mugent_ph() returns length-1 numeric components", {
  result <- mugent_ph(make_Z(5, 2), make_ldlist(2, 5))
  for (nm in names(result)) {
    expect_length(result[[nm]], 1L)
    expect_true(is.numeric(result[[nm]]))
  }
})

test_that("mugent_ph() equal Z-columns yield a large p-value (no heterogeneity)", {
  # Identical Z columns → no population heterogeneity
  set.seed(5)
  z <- rnorm(5)
  Z <- cbind(z, z)
  result <- mugent_ph(Z, make_ldlist(2, 5))
  expect_gt(result$pval, 0.05)
})

# ── mugent_pleio ──────────────────────────────────────────────────────────────

test_that("mugent_pleio() returns a named list with the expected two components", {
  result <- mugent_pleio(make_Z(5, 2), make_ldlist(2, 5))
  expect_named(result, c("result", "adjusted_significance_quantile"))
})

test_that("mugent_pleio() result field is 'pleiotropy' or 'no pleiotropy'", {
  result <- mugent_pleio(make_Z(5, 2), make_ldlist(2, 5))
  expect_true(result$result %in% c("pleiotropy", "no pleiotropy"))
})

test_that("mugent_pleio() all-zero Z-statistics return 'no pleiotropy'", {
  Z <- matrix(0, 5, 2)
  result <- mugent_pleio(Z, make_ldlist(2, 5))
  expect_equal(result$result, "no pleiotropy")
})

test_that("mugent_pleio() adjusted_significance_quantile is a positive finite scalar", {
  result <- mugent_pleio(make_Z(5, 2), make_ldlist(2, 5))
  expect_length(result$adjusted_significance_quantile, 1L)
  expect_gt(result$adjusted_significance_quantile, 0)
  expect_true(is.finite(result$adjusted_significance_quantile))
})

test_that("mugent_pleio() stricter alpha shrinks the significance quantile", {
  ldlist <- make_ldlist(2, 5)
  Z <- make_Z(5, 2)
  q_loose  <- mugent_pleio(Z, ldlist, alpha = 0.05)$adjusted_significance_quantile
  q_strict <- mugent_pleio(Z, ldlist, alpha = 1e-6)$adjusted_significance_quantile
  expect_gt(q_strict, q_loose)
})

# ── mugent_sel ────────────────────────────────────────────────────────────────

test_that("mugent_sel() returns a p x p matrix", {
  skip_if_not_installed("mvnfast")
  result <- mugent_sel(make_Z(5, 3, seed = 10), make_ldlist(3, 5, seed = 20), verbose = FALSE)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3L, 3L))
})

test_that("mugent_sel() upper triangle and diagonal are NA", {
  skip_if_not_installed("mvnfast")
  result <- mugent_sel(make_Z(5, 3, seed = 11), make_ldlist(3, 5, seed = 21), verbose = FALSE)
  expect_true(is.na(result[1, 2]))
  expect_true(is.na(result[1, 3]))
  expect_true(is.na(result[1, 1]))
  expect_true(is.na(result[2, 2]))
})

test_that("mugent_sel() lower-triangle values are in [0, 1]", {
  skip_if_not_installed("mvnfast")
  result <- mugent_sel(make_Z(5, 3, seed = 12), make_ldlist(3, 5, seed = 22), verbose = FALSE)
  lower_vals <- result[lower.tri(result)]
  lower_vals <- lower_vals[!is.na(lower_vals)]
  expect_true(all(lower_vals >= 0 & lower_vals <= 1))
})
