#Function file: montecarlo.R


# Mock .getTip function since it's not in global environment
.getTip <- function(x) x

test_that("frequencySigs identifies the most frequently occuring taxa", {
  test_sigs <- list(
    "UP_sig1" = c("Bacteroides", "Escherichia"),
    "UP_sig2" = c("Bacteroides", "Bifidobacterium")
  )
  
  result <- frequencySigs(test_sigs, n = 2)
  
  expect_s3_class(result, "table")
  expect_equal(length(result), 2)
  expect_equal(names(result)[1], "Bacteroides")
})
