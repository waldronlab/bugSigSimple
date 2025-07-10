#Function file: simple.R


# Unit test: subsetByCurator should return rows with the specified curator

library(testthat)

test_that("subsetByCurator returns only rows with the specified curator", {
  # A small fake dataset similar to what importBugSigDB() would return
  dat <- data.frame(
    Curator = c("Anne-mariesharp", "Victoria", "Anne-mariesharp"),
    Study = c("Study1", "Study2", "Study3")
  )
  
  # Use the function to subset
  result <- subsetByCurator(dat, curator = "Anne-mariesharp")
  
  # Check that the result only contains the specified curator
  expect_equal(nrow(result), 2)
  expect_true(all(result$Curator == "Anne-mariesharp"))
})
