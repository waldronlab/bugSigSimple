#Function file: simple.R

# Firstly, test when there are multiple levels
test_that(".getTip returns the last part of the string", {
  s <- "k__Bacteria|p__Firmicutes|g__Lactobacillus"
  
  result <- .getTip(s)

  expect_equal(result, "g__Lactobacillus")
})


# Secondly, test when there is only one level
test_that(".getTip works when string has only one level", {
  s <- "k__Bacteria"
  result <- .getTip(s)
  expect_equal(result, "k__Bacteria")
})
