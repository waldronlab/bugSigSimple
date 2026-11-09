# Delivery Mode and the Meconium Microbiome

## Analysis objectives

1.  Import, recode, and subset data from bugsigdb.org  
2.  Create a table of studies
3.  Calculate the frequency of appearance of each taxa in independent
    signatures and identify the most frequently reported taxa

## Making sure packages are installed

Not evaluated in vignette:

``` r

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("devtools", "tidyverse", "kableExtra"))
BiocManager::install(c("waldronlab/bugSigSimple", "waldronlab/BugSigDBStats", "waldronlab/bugsigdbr"))
```

## Load and subset data

``` r

suppressPackageStartupMessages({
  library(bugSigSimple)
  library(BugSigDBStats)
  library(bugsigdbr)
  library(tidyverse)
  library(kableExtra)
})
```

``` r

# use version="devel" and cache = FALSE to take the latest version from bugsigdb.org
dat <- bugsigdbr::importBugSigDB(cache = FALSE) 
```

``` r

dim(dat)
```

    ## [1] 7425   51

``` r

names(dat)
```

    ##  [1] "BSDB ID"                    "Study"                     
    ##  [3] "Study design"               "PMID"                      
    ##  [5] "DOI"                        "URL"                       
    ##  [7] "Authors list"               "Title"                     
    ##  [9] "Journal"                    "Year"                      
    ## [11] "Keywords"                   "Experiment"                
    ## [13] "Location of subjects"       "Host species"              
    ## [15] "Body site"                  "UBERON ID"                 
    ## [17] "Condition"                  "EFO ID"                    
    ## [19] "Group 0 name"               "Group 1 name"              
    ## [21] "Group 1 definition"         "Group 0 sample size"       
    ## [23] "Group 1 sample size"        "Antibiotics exclusion"     
    ## [25] "Sequencing type"            "16S variable region"       
    ## [27] "Sequencing platform"        "Data transformation"       
    ## [29] "Statistical test"           "Significance threshold"    
    ## [31] "MHT correction"             "LDA Score above"           
    ## [33] "Matched on"                 "Confounders controlled for"
    ## [35] "Pielou"                     "Shannon"                   
    ## [37] "Chao1"                      "Simpson"                   
    ## [39] "Inverse Simpson"            "Richness"                  
    ## [41] "Signature page name"        "Source"                    
    ## [43] "Curated date"               "Curator"                   
    ## [45] "Revision editor"            "Description"               
    ## [47] "Abundance in Group 1"       "MetaPhlAn taxon names"     
    ## [49] "NCBI Taxonomy IDs"          "State"                     
    ## [51] "Reviewer"

## Subsetting

``` r

included.pmid <-
  c(
    28018325,
    24614698,
    29207565,
    29459704,
    29538354,
    32012716,
    20566857,
    28512451,
    28112736,
    27362264
  )
subset.dat <-
  dat %>%
  dplyr::filter(PMID %in% included.pmid) 
```

``` r

unique(subset.dat$`Group 0 name`)
```

    ##  [1] "vaginal delivery"                             
    ##  [2] ">/=33 weeks"                                  
    ##  [3] "Maternal samples taken at delivery."          
    ##  [4] "Maternal samples taken at 6  weeks postpartum"
    ##  [5] "vaginal delivery at day 3"                    
    ##  [6] "vaginal delivery at day 5"                    
    ##  [7] "vaginal delivery at day 28"                   
    ##  [8] "vaginal delivery at day 150"                  
    ##  [9] "vaginal delivery at day 365"                  
    ## [10] "full-term delivery >39 weeks at day 3"        
    ## [11] "mecomium in vaginal delivery"                 
    ## [12] "mecomium in C-section delivery"               
    ## [13] "Antimicrobials use during delivery (no)"      
    ## [14] "Maternal consumption of probiotics (NO)"      
    ## [15] "Furry pets at home (NO)"

``` r

included.group0 <- "vaginal delivery"
```

``` r

unique(subset.dat$`Group 1 name`)
```

    ##  [1] "C-section"                                      
    ##  [2] "infants <33 weeks gestational age"              
    ##  [3] "Neonatal samples taken at delivery"             
    ##  [4] "Infant samples taken at 6 weeks postpartum"     
    ##  [5] "c-section"                                      
    ##  [6] "late preterm 34-36 weeks"                       
    ##  [7] "overall time point infants samples of c-section"
    ##  [8] "transitional stool"                             
    ##  [9] "Antimicrobials use during delivery (yes)"       
    ## [10] "Maternal consumption of probiotics (YES)"       
    ## [11] "Furry pets at home (YES)"

``` r

included.group1 <- "C-section"
```

``` r

subset.final <-
  subset.dat %>%
  dplyr::filter(`Group 0 name` %in% included.group0, `Group 1 name` %in% included.group1) %>%
  dplyr::filter(`Body site` == "Meconium") %>%
  dplyr::arrange(PMID)
```

Show key characteristics of the included signatures:

``` r

subset.final %>%
  dplyr::select("PMID", "Source", "Group 0 name", "Group 1 name", "Abundance in Group 1")
```

    ##       PMID              Source     Group 0 name Group 1 name
    ## 1 24614698 Table S6 in File S1 vaginal delivery    C-section
    ## 2 27362264            Figure 1 vaginal delivery    C-section
    ## 3 29459704            Figure 1 vaginal delivery    C-section
    ## 4 29459704            Figure 1 vaginal delivery    C-section
    ## 5 29538354             Table 3 vaginal delivery    C-section
    ##   Abundance in Group 1
    ## 1            increased
    ## 2            decreased
    ## 3            increased
    ## 4            decreased
    ## 5            decreased

Are any studies missing?

``` r

sort(setdiff(included.pmid, subset.dat$PMID))
```

    ## [1] 20566857 28018325 32012716

``` r

sort(setdiff(included.pmid, subset.final$PMID))
```

    ## [1] 20566857 28018325 28112736 28512451 29207565 32012716

``` r

sort(setdiff(subset.dat$PMID, subset.final$PMID))
```

    ## [1] 28112736 28512451 29207565

## Table of studies

These are the studies included in the review:

``` r

createStudyTable(subset.final) %>% 
  kbl() %>% 
  kable_styling()
```

| Study Identifier | MaxCases | MaxControls | Study design | Condition | N_signatures | PMID | DOI | URL |
|:---|---:|---:|:---|:---|---:|:---|:---|:---|
| ArdissoneA.N.\_2014 | 33 | 19 | cross-sectional observational, not case-control | Cesarean section | 1 | 24614698 | <https://doi.org/10.1371/journal.pone.0090784> | NA |
| MartinR.\_2016 | 28 | 80 | cross-sectional observational, not case-control | Cesarean section | 1 | 27362264 | <https://doi.org/10.1371/journal.pone.0158498> | NA |
| ShiY.C.\_2018 | 10 | 8 | cross-sectional observational, not case-control | Cesarean section | 2 | 29459704 | <https://doi.org/10.1038/s41598-018-21657-7> | NA |
| TapiainenT.\_2018 | 40 | 172 | cross-sectional observational, not case-control | Cesarean section | 1 | 29538354 | <https://doi.org/10.1038/pr.2018.29> | NA |

## Summary of taxa reported

This table summarizes the results for the top `n` most frequently
identified taxa.

``` r

bugSigSimple::createTaxonTable(subset.final, n = 20) %>% 
  kbl() %>% 
  kable_styling()
```

    ## Warning: Expected 7 pieces. Additional pieces discarded in 1 rows [7].

| Taxon Name | Taxonomic Level | total_signatures | increased_signatures | decreased_signatures | Binomial Test pval |
|:---|:---|---:|---:|---:|---:|
| Phocaeicola vulgatus | species | 2 | 1 | 1 | 1 |
| Actinomycetota | phylum | 1 | 0 | 1 | 1 |
| Actinomycetes | class | 1 | 0 | 1 | 1 |
| Actinomycetales | order | 1 | 0 | 1 | 1 |
| Bifidobacterium bifidum | species | 1 | 0 | 1 | 1 |
| Bifidobacterium catenulatum | species | 1 | 0 | 1 | 1 |
| Bifidobacterium longum subsp. longum | species | 1 | 0 | 1 | 1 |
| Propionibacteriaceae | family | 1 | 0 | 1 | 1 |
| Cutibacterium acnes | species | 1 | 0 | 1 | 1 |
| Propionibacterium | genus | 1 | 0 | 1 | 1 |
| Pseudonocardiaceae | family | 1 | 0 | 1 | 1 |
| Atopobium | genus | 1 | 0 | 1 | 1 |
| Bacillota | phylum | 1 | 1 | 0 | 1 |
| Bacilli | class | 1 | 1 | 0 | 1 |
| Bacillales | order | 1 | 1 | 0 | 1 |
| Bacillaceae | family | 1 | 1 | 0 | 1 |
| Bacillus | genus | 1 | 1 | 0 | 1 |
| Bacillus amyloliquefaciens | species | 1 | 1 | 0 | 1 |
| Bacillus licheniformis | species | 1 | 1 | 0 | 1 |
| Staphylococcus | genus | 1 | 0 | 1 | 1 |

## Long list of most frequently identified taxa

These are not needed because of the taxon table above, but they list a
larger number of taxa.

``` r

BugSigDBStats::getMostFrequentTaxa(subset.final, n = 50)
```

    ##    821 118964   1236   1239   1243   1279   1297   1298   1350 135614   1380 
    ##      2      1      1      1      1      1      1      1      1      1      1 
    ##   1385   1386   1390   1402   1502   1596   1598   1679   1681   1686   1730 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   1743   1747   1760 183710 186806 186817 188787 201174   2037   2070 216851 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   2737  28116  28211  28216    283    285    286  31957  32033  39491  40323 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  47678    654  72274  80840  80864    817 
    ##      1      1      1      1      1      1

``` r

BugSigDBStats::getMostFrequentTaxa(subset.final, direction="UP")
```

    ## 118964   1239   1243   1297   1298 135614   1385   1386   1390   1402 
    ##      1      1      1      1      1      1      1      1      1      1

``` r

BugSigDBStats::getMostFrequentTaxa(subset.final, direction="DOWN")
```

    ## 1236 1279 1350 1380 1502 1596 1598 1679 1681 1686 
    ##    1    1    1    1    1    1    1    1    1    1

## Acknowledgements

This vignette was refactored with the assistance of Gemini Code Assist.
