#Function file: describe_curation.R

bsdb <- bugsigdbr::importBugSigDB()
test_that("createTaxonTable creates a table of most frequent taxa in a data.frame", {
  expect_s3_class(createTaxonTable(bsdb[1:10,]), "data.frame")
})
