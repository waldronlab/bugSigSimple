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
#' @param curator curator to subset by
#' @param curator.column name of column of curators in data.frame
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
#' @param n Number of most frequently occurring taxa to return
#' @param sig.type increased for increased in cases relative to controls, decreased for decreased in cases relative to controls, both for either
#' @param direction.column column containing direction information in dat
#' @importFrom dplyr filter %>%
#' @importFrom bugsigdbr getSignatures
#' @export
#' @return 
#' A list of signatures, with PMIDs for list element names 
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
