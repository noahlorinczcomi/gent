test_that("tr() computes matrix trace correctly", {
  expect_equal(gent:::tr(diag(3)), 3)
  expect_equal(gent:::tr(matrix(c(1, 2, 3, 4), 2, 2)), 5)
  expect_equal(gent:::tr(matrix(0, 3, 3)), 0)
  expect_equal(gent:::tr(matrix(1, 1, 1)), 1)
})

test_that("f_meff() returns a positive finite value for a valid LD matrix", {
  expect_true(is.finite(gent:::f_meff(diag(5))))
  expect_gt(gent:::f_meff(diag(5)), 0)
})

test_that("f_meff() accepts both a matrix and a list of matrices", {
  ld <- diag(5)
  expect_equal(gent:::f_meff(ld), gent:::f_meff(list(ld)))
})

test_that("f_meff() returns smaller value for highly correlated SNPs vs independent", {
  m_corr <- matrix(0.99, 5, 5); diag(m_corr) <- 1
  m_indep <- diag(5)
  expect_lt(gent:::f_meff(m_corr), gent:::f_meff(m_indep))
})

test_that("f_meff() with averages two identical LD matrices equals one alone", {
  ld <- diag(5)
  expect_equal(gent:::f_meff(list(ld, ld)), gent:::f_meff(ld))
})

test_that("S() soft-thresholds positive values correctly", {
  expect_equal(gent:::S(3, 1), 2)
  expect_equal(gent:::S(1.5, 1), 0.5)
  expect_equal(gent:::S(0.5, 1), 0)
  expect_equal(gent:::S(0, 0), 0)
})

test_that("S() soft-thresholds negative values correctly", {
  expect_equal(gent:::S(-3, 1), -2)
  expect_equal(gent:::S(-0.5, 1), 0)
})

test_that("S() is vectorised", {
  expect_equal(gent:::S(c(3, -3, 0.5), 1), c(2, -2, 0))
})

test_that("pos() passes through non-negative scalars unchanged", {
  expect_equal(gent:::pos(3), 3)
  expect_equal(gent:::pos(0), 0)
})

test_that("pos() maps negative values to zero", {
  expect_equal(gent:::pos(-3), 0)
  expect_equal(gent:::pos(-1e-10), 0)
})

test_that("posadj() adjusts an indefinite matrix to be positive semi-definite", {
  # [[1, 1.2], [1.2, 1]] has eigenvalues 2.2 and -0.2 (one negative)
  r <- matrix(c(1, 1.2, 1.2, 1), 2, 2)
  expect_true(any(eigen(r)$values < 0))  # confirm it needs adjustment
  r2 <- gent:::posadj(r)
  expect_true(all(eigen(r2)$values >= -1e-10))
})

test_that("posadj() leaves an already positive-definite matrix unchanged", {
  r <- diag(3)
  expect_equal(gent:::posadj(r), r)
})

test_that("gampars() null path returns correct shape and rate for identity LD", {
  # mu = length(z) = m; sigma = 2*tr(I^2) = 2*m
  # => rate = mu/sigma = 0.5; shape = mu^2/sigma = m/2
  m <- 5
  r <- gent:::gampars(rep(0, m), diag(m), null = TRUE)
  expect_named(r, c("shape", "rate"))
  expect_equal(r$shape, m / 2)
  expect_equal(r$rate, 0.5)
})

test_that("gampars() null path is independent of z values", {
  # The null distribution parameters depend only on LD, not on z
  LD <- diag(5)
  r1 <- gent:::gampars(rep(0, 5), LD, null = TRUE)
  r2 <- gent:::gampars(rep(3, 5), LD, null = TRUE)
  expect_equal(r1$shape, r2$shape)
  expect_equal(r1$rate,  r2$rate)
})

test_that("varmat() returns a p^2 x p^2 matrix", {
  ldlist <- list(diag(5), diag(5))
  V <- gent:::varmat(2, ldlist)
  expect_equal(dim(V), c(4L, 4L))
})

test_that("varmat() stops when ldlist length does not match p", {
  expect_error(gent:::varmat(3, list(diag(5), diag(5))))
})

test_that("varmat() returns a symmetric matrix for identical LD inputs", {
  ldlist <- list(diag(5), diag(5))
  V <- gent:::varmat(2, ldlist)
  expect_equal(V, t(V))
})

test_that("varmat() diagonal entries are non-negative", {
  ldlist <- list(diag(5), diag(5))
  V <- gent:::varmat(2, ldlist)
  expect_true(all(diag(V) >= 0))
})
