#Function file: describe_curation.R

bsdb <- bugsigdbr::importBugSigDB()
test_that("createTaxonTable creates a table of most frequent taxa in a data.frame", {
  expect_s3_class(createTaxonTable(bsdb[1:10,]), "data.frame")
})


test_that("createTaxonTable reports BH-adjusted p-values", {
  tt <- createTaxonTable(bsdb[1:50, ], n = 10)
  expect_true(all(c("Binomial Test pval", "BH adjusted pval") %in% colnames(tt)))
  raw <- tt$`Binomial Test pval`
  adj <- tt$`BH adjusted pval`
  expect_true(all(adj >= raw - 1e-12, na.rm = TRUE))
  expect_true(all(adj <= 1, na.rm = TRUE))
})
