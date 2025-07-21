#Function file: simple.R


test_that("getMostFrequentTaxa returns the most frequently occurring taxa in a table of signatures", {
  # Create test data
  test_data <- data.frame(
    PMID = c("123", "456", "789"),
    `Abundance in Group 1` = c("increased", "decreased", "increased"),
    check.names = FALSE
  )
  
  # Mock bugsigdbr::getSignatures to return test taxa
  mock_signatures <- list(
    "123" = c("k__Bacteria|g__Blautia", "k__Bacteria|g__Bacteroides"),
    "456" = c("k__Bacteria|g__Blautia", "k__Bacteria|g__Lactobacillus"), 
    "789" = c("k__Bacteria|g__Blautia")
  )
  
  with_mocked_bindings(
    getSignatures = function(...) mock_signatures,
    .package = "bugsigdbr",
    {
      result <- getMostFrequentTaxa(test_data, n = 5)
      
      # Basic checks
      expect_true(is.table(result))
      expect_true(length(result) <= 5)
      expect_equal(as.numeric(result[1]), 3)  # Blautia appears 3 times
      expect_equal(names(result[1]), "k__Bacteria|g__Blautia")  # Check name
      expect_true(all(diff(result) <= 0))  # Sorted descending
    }
  )
})