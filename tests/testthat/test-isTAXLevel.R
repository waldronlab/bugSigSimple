#Function file: simple.R


test_that(".isTaxLevel returns TRUE for genus-level taxon", {
  s <- "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Lachnospira"
  result <- .isTaxLevel(s, "genus")
  expect_true(result)
})

test_that(".isTaxLevel returns FALSE for non-genus-level taxon", {
  s <- "k__Bacteria|p__Firmicutes|f__Lachnospiraceae"
  result <- .isTaxLevel(s, "genus")
  expect_false(result)
})

test_that(".isTaxLevel returns input when tax.level is 'mixed'", {
  s <- "k__Bacteria|p__Firmicutes|f__Lachnospiraceae"
  result <- .isTaxLevel(s, "mixed")
  expect_equal(result, s)
})
