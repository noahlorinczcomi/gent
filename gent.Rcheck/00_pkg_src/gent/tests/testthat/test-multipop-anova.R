make_effects <- function(m, p, seed = 1) {
  set.seed(seed)
  B  <- matrix(rnorm(m * p, sd = 0.01), m, p)
  SE <- matrix(runif(m * p, 0.005, 0.05), m, p)
  list(B = B, SE = SE)
}

test_that("multipop_anova() returns a matrix with columns ANOVA_pval and ANOVA_chisq", {
  e <- make_effects(5, 2)
  result <- multipop_anova(e$B, e$SE)
  expect_true(is.matrix(result))
  expect_equal(colnames(result), c("ANOVA_pval", "ANOVA_chisq"))
})

test_that("multipop_anova() returns one row per SNP", {
  e <- make_effects(10, 3)
  result <- multipop_anova(e$B, e$SE)
  expect_equal(nrow(result), 10L)
})

test_that("multipop_anova() p-values are in [0, 1]", {
  e <- make_effects(5, 2)
  result <- multipop_anova(e$B, e$SE)
  expect_true(all(result[, "ANOVA_pval"] >= 0))
  expect_true(all(result[, "ANOVA_pval"] <= 1))
})

test_that("multipop_anova() chi-square statistics are non-negative", {
  e <- make_effects(5, 2)
  result <- multipop_anova(e$B, e$SE)
  expect_true(all(result[, "ANOVA_chisq"] >= 0))
})

test_that("multipop_anova() stops when B and SE dimensions differ", {
  B  <- matrix(1, 5, 2)
  SE <- matrix(1, 5, 3)
  expect_error(multipop_anova(B, SE))
})

test_that("multipop_anova() identical effect columns yield large p-values", {
  # Same beta in every population → no heterogeneity → large p-value
  set.seed(1)
  b  <- rnorm(5, sd = 0.01)
  B  <- cbind(b, b)
  SE <- matrix(0.01, 5, 2)
  result <- multipop_anova(B, SE)
  expect_true(all(result[, "ANOVA_pval"] > 0.05))
})

test_that("multipop_anova() strongly opposing effects yield small p-values", {
  # Effect in one population is large and positive, large and negative in the other
  B  <- matrix(c(rep(0.5, 5), rep(-0.5, 5)), 5, 2)
  SE <- matrix(0.001, 5, 2)
  result <- multipop_anova(B, SE)
  expect_true(all(result[, "ANOVA_pval"] < 0.05))
})

test_that("multipop_anova() works with a single row (one SNP)", {
  B  <- matrix(c(0.1, -0.1), 1, 2)
  SE <- matrix(c(0.01, 0.01), 1, 2)
  result <- multipop_anova(B, SE)
  expect_equal(nrow(result), 1L)
  expect_true(is.finite(result[1, "ANOVA_pval"]))
})
