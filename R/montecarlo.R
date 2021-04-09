#' Identify the most frequently recurring taxa in bugsigdb
#'
#' @param dat data.frame produced by \link[bugsigdbr]{importBugSigDB}
#' @param direction One or both of c("increased", "decreased")
#' @param n Number of top taxa to return
#'
#' @return integer vector of top most recurrent taxa
#' @export
#'
#' @examples
#' dat <- bugsigdbr::importBugSigDB()
#' sigs <- extractBugSigs(dat, direction = "increased", tax.level = "genus")
#' frequencySigs(sigs)
frequencySigs <- function(sigs, n = 10){
  sigs.tab <- sort(table(unlist(sigs)), decreasing=TRUE)
  names(sigs.tab) <- vapply(names(sigs.tab), .getTip,
                            character(1), USE.NAMES = FALSE)
  head(sigs.tab, n=n) 
}

#' Extract a list of signatures from a bugsigdb.org table
#'
#' @param dat data.frame produced by \link[bugsigdbr]{importBugSigDB}
#' @param tax.level Any subset of c("kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")
#' @param direction One or both of c("increased", "decreased")

#' @return a list of signatures
#' @export
#'
#' @examples
#' dat <- bugsigdbr::importBugSigDB()
#' extractBugSigs(dat, tax.level = "genus")

extractBugSigs <-
  function(dat,
           direction = c("increased", "decreased"),
           tax.level = c("kingdom",
                         "phylum",
                         "class",
                         "order",
                         "family",
                         "genus",
                         "species",
                         "strain")) {
    dat <- dat[dat$`Abundance in Group 1` %in% direction, ]
    res <- extractSignatures(dat,
                             tax.level = tax.level,
                             exact.tax.level = TRUE,
                             col = "MetaPhlAn taxon names")
    zerolen <- sapply(res, function(x) length(x) == 0)
    if(any(zerolen)){
      message("Zero length signatures dropped:", paste(names(res)[zerolen], collapse = " "))
    }
    res[!zerolen]
  }

#' Simulate a list of signatures based on a universe of taxa equal in size and length to some smaller set of signatures
#'
#' @param relevant.sigs A list of signatures providing a relevant background for drawing taxa
#' @param siglengths An integer vector, the length of which provides the number of signatures to be simulated, and the 
#' integers of which provide the number of taxa to be simulated in each signature
#' @param tax.level Any subset of c("kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")
#'
#' @return A list of signatures of the same number and individual lengths as found in my.dat
#' @export
#'
#' @examples
#' full.dat <- bugsigdbr::importBugSigDB()
#' my.dat <- full.dat[full.dat$Curator == "Mst Afroza Parvin", ]
#' relevant.dat <- full.dat[full.dat$`Body site` %in% my.dat$`Body site`, ]
#' relevant.sigs <- extractBugSigs(relevant.dat, tax.level = "genus")
#' siglengths <- sapply(my.dat, length)
#' simulateSignatures(relevant.sigs, siglengths, tax.level = "genus")

simulateSignatures <-
  function(relevant.sigs,
           siglengths) {
    sigs.universe <- unlist(relevant.sigs)
    universetable <- table(sigs.universe)
    universetable <- universetable / sum(universetable)
    lapply(siglengths, function(n)
      sample(names(universetable), size = n, prob = universetable))
  }

.countBug <- function(relevant.sigs, siglengths){
  siglist <- simulateSignatures(relevant.sigs, siglengths)
  max(table(unlist(siglist)))
}

#' countBug counts the frequency of the most commonly identified bug in a simulated signature.

#' getCriticalN performs a Monte Carlo simulation to estimate the number of times the most frequent taxon is expected to be observed
#' in a list of signatures
#'
#' @param relevant.sigs a list of signatures that form the "background" from which taxa for simulated signatures will be drawn. 
#' These are used to estimate how frequently taxa occur
#' @param siglengths The sizes of signatures found in a set of related studies. Simulated signatures will match these in number and size.
#' @param alpha Probability at which a critical threshold will be calculated (default: 0.05)
#' @param nsim Number of simulations (default: 1000)
#'
#' @return The 1 - alpha quantile of Monte Carlo simulated values for the maximum number of times any taxon is identified.
#' @export
#' @details E.g. for alpha = 0.05, we expect only a 5% chance that any taxon will be identified N times or more.

#' @examples
#' full.dat <- bugsigdbr::importBugSigDB()
#' my.dat <- full.dat[full.dat$Curator == "Mst Afroza Parvin", ]
#' relevant.dat <- full.dat[full.dat$`Body site` %in% my.dat$`Body site`, ]
#' relevant.sigs <- extractBugSigs(my.dat)
#' my.sigs.increased <- extractBugSigs(my.dat, direction = "increased")
#' (my.siglengths <- sapply(my.sigs.increased, length))
#' getCriticalN(relevant.sigs, my.siglengths)
#' # Compare to observed
#' frequencySigs(my.sigs.increased)

getCriticalN <- function(relevant.sigs, siglengths, alpha = 0.05, nsim = 1000){
  res <- replicate(nsim, suppressMessages(.countBug(relevant.sigs, siglengths)))
  quantile(res, 1 - alpha)
}