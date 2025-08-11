#' Create a table of most frequent taxa in a data.frame
#'
#' @param dat data.frame produced by \link[bugsigdbr]{importBugSigDB}, subsetted as desired
#' @param n number of taxa to return (if sig.type=="both", this is the number of taxa to return for each direction)
#'
#' @importFrom dplyr filter %>% mutate rowwise n_distinct group_by relocate rename first ungroup across n
#' @importFrom utils relist head
#' @importFrom stats quantile binom.test
#' @importFrom kableExtra kbl kable_styling
#' @importFrom tidyr separate
#' @importFrom stringr str_replace str_extract
#' @importFrom bugsigdbr getSignatures
#' 
#' @return kable table with increased and decreased taxa and a binomial test based on total number of studies in the data.frame
#' @export
#'
#' @examples
#' full.dat <- bugsigdbr::importBugSigDB()
#' createTaxonTable(full.dat, n=20)

createTaxonTable <- function(dat, n=10){
  dmap <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  names(dmap) <- substring(dmap, 1, 1)
  output <-
    data.frame(getMostFrequentTaxa(dat, sig.type = "both", n = n),
               stringsAsFactors = FALSE) %>%
    mutate(metaphlan_name = Var1) %>%
    separate(
      col = Var1,
      sep = "\\|",
      into = dmap,
      fill = "right"
    ) %>%
    mutate(across(kingdom:species, ~ str_replace(., ".__", ""))) %>%
    rename(n_signatures = Freq)
  
  output <-
    output %>% mutate(n_signatures = sapply(output$metaphlan_name, function(x) {
      sum(grepl(
        pattern = x,
        x = dat$`MetaPhlAn taxon names`,
        fixed = TRUE
      ))
    })) %>%
    mutate(`Total Signatures` = sapply(metaphlan_name, function(x)
      .countTaxon(dat = dat, x = x, direction = "both"))) %>%
    mutate(`Increased Signatures` = sapply(metaphlan_name, function(x)
      .countTaxon(dat = dat, x = x, direction = "increased"))) %>%
    mutate(`Decreased Signatures` = sapply(metaphlan_name, function(x)
      .countTaxon(dat = dat, x = x, direction = "decreased"))) %>%
    mutate(Taxon = gsub(".+\\|", "", output$metaphlan_name))
  
    output %>% separate(col="Taxon", into=c("Taxonomic Level", "Taxon Name"), sep="__") %>%
    mutate(`Taxonomic Level` = unname(dmap[`Taxonomic Level`])) %>%
    rowwise() %>%    
    mutate( `Binomial Test pval` = .createBinomTestSummary(`Increased Signatures`, `Total Signatures`, wordy = FALSE)) %>%
      ungroup() %>%
    relocate(`Taxon Name`, `Taxonomic Level`, `Total Signatures`, `Increased Signatures`, `Decreased Signatures`, `Binomial Test pval`)
}

.countTaxon = function(dat, x, direction = c("both", "increased", "decreased")){
  if (direction[1] %in% c("increased", "decreased")){
    dat <- filter(dat, `Abundance in Group 1` == direction[1])
  }
  allnames <- getSignatures(dat, tax.id.type = "metaphlan")
  sum(vapply(allnames, function(onesignames) x %in% onesignames, FUN.VALUE = 1L))
}

.createBinomTestSummary <- function(x, n, p=0.5, wordy = TRUE){
  bin.test <- binom.test(x=x, n=n, p=p)
  if(wordy){
    paste0("p-value=", signif(bin.test[[3]],2), " (increased freq=", bin.test[[1]],", total freq=", bin.test[[2]], ")") %>% return()
  }else{
    signif(bin.test[[3]], 2)
  }
}


.getLowestTaxon <- function(taxa){
  taxaNames <- taxa %>% names() 
  tax.level <- lapply(taxaNames, function(x)
  {
    spl <- unlist(strsplit(x, "\\|"))
    leave <- spl[length(spl)]
  })
  tax.level %>% unlist() %>% data.frame(`Taxon`=., taxa) %>% return()
}

#' Create a table of all studies currently in data.frame
#'
#' @param bsdb.df \code{data.frame} produced by \link[bugsigdbr]{importBugSigDB}, pre-filtered as desired
#' @param includeAlso \code{character} with column names to additionally include in the output table (default = `NULL`)
#'
#' @importFrom dplyr group_by %>% select relocate reframe across n
#' @importFrom tidyr all_of
#' @importFrom kableExtra kbl kable_styling
#' @return a data.frame of basic study information. Can be wrapped in 
#' kable_styling(kbl(...)) to format nicely.
#' @export
#'
#' @examples
#' full.dat <- bugsigdbr::importBugSigDB()
#' createStudyTable(full.dat)

createStudyTable <- function(bsdb.df, includeAlso = NULL) {
  # input check
  if (!is.null(includeAlso)) {
    if (!all(includeAlso %in% colnames(bsdb.df))) {
      stop(paste(
        "The following columns are not found in the input data frame:",
        paste(includeAlso[!(includeAlso %in% colnames(bsdb.df))], collapse = ", ")
      ))
    }
  }
  # Core of the change is in how study IDs are generated, see function in 
  # simple.R. NB: the function also fixes DOI links as side effect, now. 
  
  bsdb_with_StudyIDs.df <- .make_unique_study_ID(bsdb.df)
  
  # some dplyr-fu to summarize tables, with more recent syntax
  study_table_fixed <- bsdb_with_StudyIDs.df %>%
    group_by(`Study Identifier`) %>%
    reframe(
      Cases = max(`Group 1 sample size`),
      Controls = max(`Group 0 sample size`),
      across(
        all_of(
          c("Study design", "Condition", "PMID", "DOI", "URL", includeAlso)
        ),
        .fns = function(x)
          paste(unique(x), collapse = "; ")
      ),
      `Number of signatures` = n()
    ) %>%
    relocate(`Number of signatures`, .after = Condition)
  
  return(study_table_fixed)
}

globalVariables(
  c(
    ".",
    "Study",
    "Condition",
    "Cases",
    "Controls",
    "Study.Design",
    "Taxon Name",
    "Binomial Test pval",
    "Total Signatures",
    "Abundance in Group 1",
    "Decreased Signatures",
    "Increased Signatures",
    "Taxonomic Level",
    "metaphlan_name",
    "Freq",
    "species",
    "kingdom",
    "Var1",
    # the following are necessary after the merged data
    "Authors",
    "Year",
    "BasicID",
    "PMID",
    "URL",
    "uniqueRank",
    "Study Identifier",
    "Group 0 sample size",
    "Group 1 sample size",
    "Number of signatures"
  )
)

