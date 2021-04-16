#' Create a table of most frequent taxa in a data.frame
#'
#' @param dat data.frame produced by \link[bugsigdbr]{importBugSigDB}, subsetted as desired
#' @param sig.type increased for increased in cases relative to controls, decreased for decreased in cases relative to controls, both for either
#' @param n number of taxa to return (if sig.type=="both", this is the number of taxa to return for each direction)
#'
#' @importFrom dplyr filter %>%
#' @importFrom kableExtra kbl kable_styling
#' @importFrom tidyr separate
#' @return kable table
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
  
  taxa <- separate(data=taxa, col="Taxon", into=c("Taxonomic Level", "Taxon Name"), sep="__")
  dmap <- c("class", "phylum", "order", "family", "species", "genus")
  names(dmap) <- substring(dmap, 1, 1)
  taxa$`Taxonomic Level` <- unname(dmap[taxa$`Taxonomic Level`])
  
  taxa %>% kbl() %>% kable_styling()
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
