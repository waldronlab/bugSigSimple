############################################################
# 
# author: Ludwig Geistlinger
# date: 2019-03-13 17:01:25
# 
# descr: simple exploratory analysis for curated signatures
# 
###########################################################

readCurationSheet <- function(data.file)
{
    dat <- read.delim(data.file, skip=1, as.is=TRUE)
    ind <- dat[,"PMID"] != ""
    dat <- dat[ind,]
    dat <- dat[-1,]
    return(dat)
}

subsetByCondition <- function(dat, condition, condition.column="condition")
{
    ind <- dat[,condition.column] == condition
    return(dat[ind,])
}

subsetByCurator <- function(dat, curator, curator.column="curator")
{
    ind <- dat[,curator.column] == curator
    return(dat[ind,])
}

extractSignatures <- function(dat)
{
    ind <- grep("^taxon.1$", colnames(dat))
    ind <- ind:ncol(dat)
    msc <- apply(as.matrix(dat[,ind]), 1, function(x) x[!is.na(x) & x != ""])
    return(msc)
}

getMostFrequentTaxa <- function(dat, n=10, sig.type=c("BOTH", "UP", "DOWN"))
{
    sig.type <- match.arg(sig.type)
    if(sig.type %in% c("UP", "DOWN")) 
    {
        ind <- dat[,"UP.or.DOWN"] == sig.type
        dat <- dat[ind,]
    }
    msc <- extractSignatures(dat)
    msc.tab <- sort(table(unlist(msc)), decreasing=TRUE)
    head(msc.tab, n=n) 
}
