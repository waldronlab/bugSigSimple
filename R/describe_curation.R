#' Create a table of most frequent taxa in a data.frame
#'
#' @param dat data.frame produced by \link[bugsigdbr]{importBugSigDB}, subsetted as desired
#' @param sig.type increased for increased in cases relative to controls, decreased for decreased in cases relative to controls, both for either
#' @param n number of taxa to return (if sig.type=="both", this is the number of taxa to return for each direction)
#'
#' @importFrom dplyr filter %>%
#' @importFrom kableExtra kbl kable_styling
#' @importFrom tidyr separate
#' @return kable table with increased and decreased taxa and a binomial test based on total number of studies in the data.frame
#' @export
#'
#' @examples
#' full.dat <- bugsigdbr::importBugSigDB()
#' createTaxonTable(full.dat, n=20)


createTaxonTable <- function(dat, sig.type=c("increased, decreased, both"), n=10){
  inc <- dat %>% getMostFrequentTaxa(., sig.type = "increased", n=n) %>% .getLowestTaxon(.) %>% {.[,-2]}
  inc <- cbind(direction=rep("increased", nrow(inc)),inc)
  dec <- dat %>% getMostFrequentTaxa(., sig.type = "decreased", n=n) %>% .getLowestTaxon(.) %>% {.[,-2]}
  dec <- cbind(direction=rep("decreased", nrow(dec)),dec)

  taxa <- rbind(inc,dec)
  if(sig.type %in% c("increased", "decreased")){taxa <- filter(taxa,direction==sig.type)}
  
  taxa <- tidyr::separate(data=taxa, col="Taxon", into=c("Taxonomic Level", "Taxon Name"), sep="__")
  dmap <- c("class", "phylum", "order", "family", "species", "genus")
  names(dmap) <- substring(dmap, 1, 1)
  taxa$`Taxonomic Level` <- unname(dmap[taxa$`Taxonomic Level`])
  nstudy <- dat$PMID %>% n_distinct()
  taxa <- taxa  %>% rowwise() %>% mutate(`Binomial Test`=.createBinomTestSummary(Freq, nstudy))
  taxa %>% kbl() %>% kable_styling()
}

.createBinomTestSummary <- function(x, n, p=0.5){
  bin.test <- binom.test(x=x, n=n, p=p)
  paste0("p-value=", round(bin.test[[3]],4), " (x=", bin.test[[1]],", n=", bin.test[[2]], ")") %>% return()
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
#' @param dat data.frame produced by \link[bugsigdbr]{importBugSigDB}, subsetted as desired
#'
#' @importFrom dplyr group_by summarize  %>%
#' @importFrom kableExtra kbl kable_styling
#' @return kable table of basic study information
#' @export
#'
#' @examples
#' full.dat <- bugsigdbr::importBugSigDB()
#' createStudyTable(full.dat)

createStudyTable <-function(dat){
  studies <- data.frame(Study=paste0(str_extract(dat$Authors, "[A-Za-z]+[:space:]"), dat$Year),
                        Condition=dat$Condition,
                        Cases=dat$`Group 1 sample size`,
                        Controls=dat$`Group 0 sample size`,
                        `Study Design`=dat$`Study design`)
  studies %>% group_by(Study) %>% summarize(Condition=first(Condition), 
                                            Cases=max(Cases),
                                            Controls=max(Controls), 
                                            `Study Design`=first(`Study.Design`)) %>%
    kbl() %>% 
    kable_styling()
  
}
