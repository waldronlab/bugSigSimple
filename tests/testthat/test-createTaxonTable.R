# ===== Load dependencies =====
library(dplyr)
library(tidyr)
library(stringr)
library(testthat)

# ===== Mock getSignatures function =====
getSignatures <- function(dat, tax.id.type = "metaphlan") {
  strsplit(dat$`MetaPhlAn taxon names`, "\\|")
}

# ===== Helper: countTaxon =====
.countTaxon <- function(dat, x, direction = c("both", "increased", "decreased")) {
  direction <- match.arg(direction)
  if (direction != "both") {
    dat <- dplyr::filter(dat, `Abundance in Group 1` == direction)
  }
  allnames <- getSignatures(dat)
  sum(vapply(allnames, function(sig) x %in% sig, integer(1)))
}

# ===== Helper: Binomial Test Summary =====
.createBinomTestSummary <- function(x, n, p = 0.5, wordy = TRUE) {
  bt <- binom.test(x, n, p)
  if (wordy) {
    paste0("p-value=", signif(bt$p.value, 2), 
           " (increased freq=", x, ", total freq=", n, ")")
  } else {
    signif(bt$p.value, 2)
  }
}

# ===== Mock getMostFrequentTaxa function =====
getMostFrequentTaxa <- function(dat, sig.type = "both", n = 10) {
  tab <- table(dat$`MetaPhlAn taxon names`)
  df <- as.data.frame(head(sort(tab, decreasing = TRUE), n))
  names(df) <- c("Var1", "Freq")
  df
}

# ===== Main Function: createTaxonTable =====
createTaxonTable <- function(dat, n = 10) {
  dmap <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  names(dmap) <- substring(dmap, 1, 1)
  
  output <- data.frame(getMostFrequentTaxa(dat, sig.type = "both", n = n),
                       stringsAsFactors = FALSE) %>%
    mutate(metaphlan_name = Var1) %>%
    separate(col = Var1, sep = "\\|", into = dmap, fill = "right") %>%
    mutate(across(kingdom:species, ~ str_replace(., ".__", ""))) %>%
    rename(n_signatures = Freq)
  
  output <- output %>%
    mutate(n_signatures = sapply(metaphlan_name, function(x) {
      sum(grepl(x, dat$`MetaPhlAn taxon names`, fixed = TRUE))
    })) %>%
    mutate(total_signatures = sapply(metaphlan_name, function(x)
      .countTaxon(dat, x, "both"))) %>%
    mutate(increased_signatures = sapply(metaphlan_name, function(x)
      .countTaxon(dat, x, "increased"))) %>%
    mutate(decreased_signatures = sapply(metaphlan_name, function(x)
      .countTaxon(dat, x, "decreased"))) %>%
    mutate(Taxon = gsub(".+\\|", "", metaphlan_name))
  
  output %>%
    separate(Taxon, into = c("Taxonomic Level", "Taxon Name"), sep = "__") %>%
    mutate(`Taxonomic Level` = unname(dmap[`Taxonomic Level`])) %>%
    rowwise() %>%
    mutate(`Binomial Test pval` = .createBinomTestSummary(increased_signatures, total_signatures, wordy = FALSE)) %>%
    ungroup() %>%
    relocate(`Taxon Name`, `Taxonomic Level`, total_signatures,
             increased_signatures, decreased_signatures, `Binomial Test pval`)
}

# ======= TEST =========

test_that("createTaxonTable returns correct structure", {
  mock_data <- data.frame(
    `MetaPhlAn taxon names` = c(
      "k__Bacteria|p__Firmicutes|g__Lactobacillus",
      "k__Bacteria|p__Firmicutes|g__Lactobacillus",
      "k__Bacteria|p__Bacteroidetes|g__Bacteroides"
    ),
    `Abundance in Group 1` = c("increased", "decreased", "increased"),
    stringsAsFactors = FALSE
  )
  
  result <- createTaxonTable(mock_data, n = 2)
  
  expect_s3_class(result, "data.frame")
  expect_true(all(c("Taxon Name", "Taxonomic Level", "total_signatures",
                    "increased_signatures", "decreased_signatures", "Binomial Test pval") %in% colnames(result)))
  expect_gt(nrow(result), 0)
})
