#Function file: simple.R


# Unit test: MPA.TAX.LEVELS should return the correct taxonomic level abbreviation

library(testthat)

test_that("MPA.TAX.LEVELS contains correct taxonomic level abbreviations", {
  expected <- c("k", "p", "c", "o", "f", "g", "s", "t")
  names(expected) <- c("kingdom", "phylum", "class", "order",
                       "family", "genus", "species", "strain")
  
  expect_equal(MPA.TAX.LEVELS, expected)
})

