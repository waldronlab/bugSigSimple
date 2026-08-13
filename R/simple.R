############################################################
# 
# author: Ludwig Geistlinger
# date: 2019-03-13 17:01:25
# 
# descr: simple exploratory analysis for curated signatures
# 
###########################################################

TAX.LEVELS <- c("kingdom", "phylum", "class", "order",
                "family", "genus", "species", "strain")
MPA.TAX.LEVELS <- c(substring(TAX.LEVELS[1:7], 1, 1), "t")
names(MPA.TAX.LEVELS) <- TAX.LEVELS

MPA.REGEXP <- "^[kpcofgst]__"

#' Subset a data.frame of signatures by curator
#'
#' @param dat data.frame produced by \link[bugsigdbr]{importBugSigDB}
#' @param curator character vector of one or more curator names to subset by
#' @param curator.column name of column of curators in data.frame (default: "Curator")
#'
#' @return data.frame subsetted by curator
#' @importFrom dplyr filter %>%
#' @export
#'
#' @examples
#' full.dat <- bugsigdbr::importBugSigDB()
#' fatima.dat <- subsetByCurator(full.dat, curator="Fatima Zohra")

subsetByCurator <- function(dat, curator, curator.column="Curator")
{
    dat %>% filter(!!as.name(curator.column) %in% !!curator) %>% return()
}

.isTaxLevel <- function(s, tax.level)
{
    if(tax.level[1] == "mixed") return(s)
    tip <- .getTip(s)
    tip <- substring(tip, 1, 1)
    mtl <- MPA.TAX.LEVELS[tax.level]
    tip %in% mtl
}

.getTip <- function(n)
{
    spl <- unlist(strsplit(n, "\\|"))
    spl[length(spl)]
}


#' Get the most frequently occurring taxa in a table of signatures
#' @param dat A table such as output by \link[bugsigdbr]{importBugSigDB}
#' @param n Number of most frequently occurring taxa to return (default: 10)
#' @param sig.type "increased" for taxa increased in cases relative to controls, "decreased" for decreased, "both" for either (default: "both")
#' @param direction.column column containing direction information in dat (default: "Abundance in Group 1")
#' @importFrom dplyr filter %>%
#' @importFrom bugsigdbr getSignatures
#' @export
#' @return
#' A named table of taxon counts (as produced by \code{table}), sorted in
#' decreasing order. Names are the full metaphlan-style taxon path strings
#' (e.g. "k__Bacteria|g__Blautia", not just the tip taxon), and values are
#' their frequency across signatures in \code{dat}.
#' @examples
#' full.dat <- bugsigdbr::importBugSigDB()
#' getMostFrequentTaxa(full.dat)
getMostFrequentTaxa <- function(dat, n=10, sig.type=c("both", "increased", "decreased"), direction.column="Abundance in Group 1")
{
    sig.type <- match.arg(sig.type)
    
    if(sig.type %in% c("increased", "decreased")) 
    {
        dat <- dat %>% filter(!!as.name(direction.column) == !!sig.type)
    }
    msc <- bugsigdbr::getSignatures(dat, tax.id.type = "metaphlan")
    msc.tab <- sort(table(unlist(msc)), decreasing=TRUE)
    head(msc.tab, n=n) 
}


#' Build unique study identifiers for a table of signatures
#'
#' Generates a study identifier the way curatedMetagenomicData does it: the
#' first author from the \code{Authors list} column (usually surname +
#' initials, e.g. "ZellerG") with whitespace stripped, plus the year, e.g.
#' "ZellerG_2014". Duplicates get a ".1", ".2", ... suffix. Also normalizes
#' bare DOIs (starting with "10.") to full \url{https://doi.org/} links.
#'
#' @param bsdb.df \code{data.frame} produced by \link[bugsigdbr]{importBugSigDB}, pre-filtered as desired
#'
#' @return \code{bsdb.df} with a new leading \code{Study Identifier} column
#' @note Author: Giacomo Antonello (2025-03-17)
#' @keywords internal
#'
#' @examples
#' full.dat <- bugsigdbr::importBugSigDB()
#' bugSigSimple:::.make_unique_study_ID(full.dat[1:10, ])

.make_unique_study_ID <- function(bsdb.df){
  bsdb_with_StudyIDs <- bsdb.df %>% 
    # fix DOIs
    mutate(
      DOI =  ifelse(
        test = startsWith(DOI, "10."),
        yes = paste0("https://doi.org/", DOI),
        no = DOI
      ),
      # create a basic ID
      BasicID = paste0(gsub(" ", "", sapply(strsplit(`Authors list`, ", "), "[", 1)), "_", Year)
    ) %>% 
    # For each ID found, seach if there are multiple studies
    group_by(BasicID) %>% 
    mutate(
      # this is arbitrary, the point is to make sure you can split overlapping
      # IDs into one
      uniqueRank = as.numeric(as.factor(paste(PMID, DOI, URL, `Authors list`))),
      `Study Identifier` = ifelse(uniqueRank > 1, paste(BasicID, uniqueRank - 1, sep = "."), BasicID)
    ) %>% 
    ungroup() %>% 
    select(- BasicID, - uniqueRank) %>% 
    relocate(`Study Identifier`) 
  
  return(bsdb_with_StudyIDs)
}
