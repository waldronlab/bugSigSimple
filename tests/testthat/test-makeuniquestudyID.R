#Function file: simple.R


test_that("make_unique_study_ID works", {
  # Test data
  test_df <- data.frame(
    DOI = c("10.1234/example", "10.5678/duplicate"),
    `Authors list` = c("Smith, John", "Smith, John"),
    Year = c(2020, 2020),
    PMID = c("111", "222"),
    URL = c("url1", "url2"),
    check.names = FALSE
  )
  
  result <- .make_unique_study_ID(test_df)
  
  # Basic checks
  expect_true("Study code" %in% names(result))
  expect_equal(nrow(result), 2)
  expect_true(all(startsWith(result$DOI, "https://doi.org/")))
  expect_length(unique(result$`Study code`), 2) # Should create unique IDs
})