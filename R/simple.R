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


#' Subset a data.frame of signatures by condition of interest
#'
#' @param dat data.frame produced by \link[bugsigdbr]{importBugSigDB}
#' @param condition health condition or disease of interest to subset by
#' @param condition.column name of column of conditions in data.frame
#'
#' @importFrom dplyr filter %>%
#' @return data.frame subsetted by condition
#' @export
#'
#' @examples
#' full.dat <- bugsigdbr::importBugSigDB()
#' obese.dat <- subsetByCondition(full.dat, condition="obesity")

subsetByCondition <- function(dat, condition, condition.column="Condition")
{
    dat %>% filter(!!as.name(condition.column) %in% condition) %>% return()
}

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

#' Create a list of signatures
#' @param dat A table such as output by \link[bugsigdbr]{importBugSigDB}
#' 
#' @param tax.level Either "mixed" or any subset of c("kingdom", "phylum", "class", "order", "family", "genus", "species", "strain"). This full vector is equivalent to "mixed". 
#' @param exact.tax.level If TRUE, return only the exact taxonomic levels specified by tax.level. FALSE is not working.
#' @param col Column name containing signatures. "MetaPhlAn taxon names" for bugsigdb.org. 
#'
#' @export
#' @return 
#' A list of signatures, with PMIDs for list element names 
#' @examples 
#' full.dat <- bugsigdbr::importBugSigDB()
#' extractSignatures(full.dat, tax.level="genus", exact.tax.level=TRUE, col = "MetaPhlAn taxon names")

#Will eventually be imported from bugSigDB package but it's not there yet
extractSignatures <- function(dat, tax.level = "mixed", 
                              exact.tax.level = TRUE, col = "MetaPhlAn taxon names")
{
    stopifnot(is.character(tax.level))
    if("mixed" %in% tax.level) tax.level <- "mixed"
    else if(!all(tax.level %in% TAX.LEVELS))
        stop("tax.level must be a subset of { ", 
             paste(TAX.LEVELS, collapse = ", "), " }")
    
    if(is.null(col))
    {    
        ind <- grep("^taxon.1$", colnames(dat))
        ind <- ind:ncol(dat)
        msc <- apply(as.matrix(dat[,ind]), 1, function(x) x[!is.na(x) & x != ""])
    }
    else msc <- strsplit(dat[,col], ",")    
    
    msc <- lapply(msc, unique)
    n <- rle(dat[,"PMID"])$lengths
    n <- unlist(lapply(n, seq_len))
    names(msc) <- paste0("PMID:", dat[,"PMID"], "_", n)
    
    dupl.ind <- duplicated(names(msc))
    dupl.pmids <- names(msc)[dupl.ind]
    dupl.pmids <- sub("_[0-9]+$", "", dupl.pmids) 
    dupl.pmids <- sub("^PMID:", "", dupl.pmids) 
    dupl.pmids <- unique(dupl.pmids)
    for(pmid in dupl.pmids)
    {
        pstr <- paste0("PMID:", pmid, "_")
        ind <- grep(pstr, names(msc))
        names(msc)[ind] <- paste0(pstr, seq_along(ind))
    } 
    
    if(tax.level[1] != "mixed")
    {
        if(!exact.tax.level)
            msc <- lapply(msc, .extractTaxLevelSig, tax.level = tax.level)
        
        bugs <- unique(unlist(msc))
        ind <- lapply(msc, function(s) match(s, bugs))
        istl <- vapply(bugs, .isTaxLevel, logical(1), tax.level = tax.level)
        subind <- istl[unlist(ind)]
        subind <- relist(subind, ind)
        msc <- mapply(`[`, msc, subind)
    }
    return(msc) 
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
    msc <- extractSignatures(dat)
    msc.tab <- sort(table(unlist(msc)), decreasing=TRUE)
    head(msc.tab, n=n) 
}
