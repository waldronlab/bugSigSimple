#Function file: simple.R


test_that("subsetByCurator returns only rows with the specified curator", {
  dat <- data.frame(
    Curator = c("Anne-mariesharp", "Victoria", "Anne-mariesharp"),
    Study = c("Study 1", "Study 2", "Study 3")
  )
  
  # Use the function to subset
  result <- subsetByCurator(dat, curator = "Anne-mariesharp")
  
  # Check that the result only contains the specified curator
  expect_equal(nrow(result), 2)
  expect_true(all(result$Curator == "Anne-mariesharp"))
})