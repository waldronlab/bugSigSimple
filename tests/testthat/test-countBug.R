#Function file: montecarlo.R


test_that(".countBug returns max frequency of most common taxon", {
  relevant.sigs <- list(
    sig1 = c("Bacteroides", "Escherichia"),
    sig2 = c("Bacteroides", "Bifidobacterium")
  )
  
  set.seed(123)
  result <- .countBug(relevant.sigs, c(2, 2))
  
  expect_type(result, "integer")
  expect_true(result >= 1)
  expect_true(result <= 4)  # Max possible count
})