#Function file: montecarlo.R


library(testthat)

test_that("simulateSignatures simulates signatures by sampling taxa based on 
          their frequency in a reference set", {
  relevant.sigs <- list(
    sig1 = c("Bacteroides", "Escherichia"),
    sig2 = c("Bacteroides", "Bifidobacterium")
  )
  
  result <- simulateSignatures(relevant.sigs, c(2, 1))
  
  expect_type(result, "list")
  expect_equal(length(result), 2)
  expect_equal(length(result[[1]]), 2)
  expect_equal(length(result[[2]]), 1)
})