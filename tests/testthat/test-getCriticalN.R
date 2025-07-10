#Function file: montecarlo.R


library(testthat)

test_that("getCriticalN returns 95th percentile threshold", {
  relevant.sigs <- list(
    sig1 = c("Bacteroides", "Escherichia"),
    sig2 = c("Bacteroides", "Bifidobacterium")
  )
  
  set.seed(123)
  result <- getCriticalN(relevant.sigs, c(2, 2), nsim = 100)
  
  expect_type(result, "double")
  expect_true(result >= 1)
  expect_true(result <= 4)
})