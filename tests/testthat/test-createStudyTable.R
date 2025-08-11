#Function file: describe_curation.R

bsdb <- bugsigdbr::importBugSigDB()
test_that("createStudyTable creates a table of all studies currently in data.frame", {
  expect_s3_class(createStudyTable(bsdb[1:10,]), "data.frame")
})