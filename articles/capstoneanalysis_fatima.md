# An analysis of Most frequent taxa in Major Depression and Bipolar Disorder

## Setup phase: Install and load relevant packages

``` r

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::version()
BiocManager::install("remotes", dependencies = TRUE)
BiocManager::install("waldronlab/bugSigSimple")
```

Load bugSigSimple and bugSigDB

``` r

library(bugSigSimple)
library(bugsigdbr)
library(kableExtra)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.2.1     ✔ readr     2.2.0
    ## ✔ forcats   1.0.1     ✔ stringr   1.6.0
    ## ✔ ggplot2   4.0.3     ✔ tibble    3.3.1
    ## ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
    ## ✔ purrr     1.2.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter()     masks stats::filter()
    ## ✖ dplyr::group_rows() masks kableExtra::group_rows()
    ## ✖ dplyr::lag()        masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

## Load latest BugSigDB data

``` r

dat <- importBugSigDB()
```

    ## Using cached version from 2026-11-09 16:29:54

``` r

dim(dat)
```

    ## [1] 7425   51

## Subset all signatures by curator and condition

``` r

curatorName <- "Fatima Zohra"

caption <- paste0("Overview of the studies initially selected based on curator name: ", curatorName, ". For convenience, only top 20 most frequently reported conditions are shown")
```

``` r

my.dat <- subsetByCurator(dat, curator = curatorName)

table(my.dat[,"Condition"]) %>% 
  sort(decreasing = TRUE) %>% 
  as.data.frame() %>% 
  head(20) %>% 
  kbl(caption = caption) %>% 
  kable_styling()
```

| Var1                                     | Freq |
|:-----------------------------------------|-----:|
| Schizophrenia                            |   19 |
| HIV infection                            |   17 |
| Bipolar disorder                         |   16 |
| Lung cancer                              |   16 |
| Gastric cancer                           |   12 |
| Chronic kidney disease                   |    8 |
| Pancreatic carcinoma                     |    8 |
| Unipolar depression                      |    8 |
| Anorexia nervosa                         |    6 |
| Colorectal cancer                        |    6 |
| Polycystic ovary syndrome                |    6 |
| Anxiety disorder                         |    5 |
| Autism                                   |    5 |
| Attention deficit hyperactivity disorder |    4 |
| Cervical cancer                          |    4 |
| Chronic hepatitis C virus infection      |    4 |
| Gastric adenocarcinoma                   |    4 |
| Hepatocellular carcinoma                 |    4 |
| Human papilloma virus infection          |    4 |
| Multiple sclerosis                       |    4 |

Overview of the studies initially selected based on curator name: Fatima
Zohra. For convenience, only top 20 most frequently reported conditions
are shown {.table .table style="margin-left: auto; margin-right: auto;"}

``` r

efo <- bugsigdbr::getOntology("efo")
```

    ## Loading required namespace: ontologyIndex

    ## Using cached version from 2026-11-09 16:21:56

``` r

dat.bpd <- bugsigdbr::subsetByOntology(dat, column = "Condition", "bipolar disorder", efo)
dat.upd <- bugsigdbr::subsetByOntology(dat, column = "Condition", "major depressive disorder", efo)
my.dat.cond <- rbind(dat.bpd, dat.upd)
table(my.dat.cond[,"Condition"])
```

    ## 
    ##          Bipolar disorder Major depressive disorder 
    ##                        18                        32

## Tables for selected conditions together

### Studies table

``` r

createStudyTable(my.dat.cond) %>% 
  kbl() %>% 
  kable_styling()
```

| Study Identifier | MaxCases | MaxControls | Study design | Condition | N_signatures | PMID | DOI | URL |
|:---|---:|---:|:---|:---|---:|:---|:---|:---|
| ChenJ.J.\_2020 | 45 | 44 | cross-sectional observational, not case-control | Major depressive disorder | 7 | 32040443 | <https://doi.org/10.18632/aging.102775> | <https://pubmed.ncbi.nlm.nih.gov/32040443/> |
| ChenY.\_2024 | 29 | 33 | case-control | Major depressive disorder | 4 | 38436041 | <https://doi.org/10.2147/NDT.S449224> | NA |
| CoelloK.\_2019 | 113 | 39 | cross-sectional observational, not case-control | Bipolar disorder | 1 | 30261302 | <https://doi.org/10.1016/j.bbi.2018.09.026> | NA |
| EvansS.J.\_2017 | 115 | 64 | case-control | Bipolar disorder | 1 | 27988330 | <https://doi.org/10.1016/j.jpsychires.2016.12.007> | NA |
| FlowersS.A.\_2017 | 49 | 68 | cross-sectional observational, not case-control | Bipolar disorder | 2 | 28035686 | <https://doi.org/10.1002/phar.1890> | NA |
| JiangY.\_2024 | 65 | 61 | prospective cohort | Major depressive disorder | 2 | 39695082 | <https://doi.org/10.1038/s41398-024-03208-z> | <https://www.nature.com/articles/s41398-024-03208-z> |
| LiY.\_2024 | NA | NA | laboratory experiment | Major depressive disorder | 13 | 38344423 | <https://doi.org/10.2147/NDT.S448940> | NA |
| LiuP.\_2024 | 85 | 85 | case-control,laboratory experiment | Major depressive disorder | 5 | 38589368 | <https://doi.org/10.1038/s41467-024-47273-w> | NA |
| Malan-MullerS.\_2022 | 20 | 2 | case-control | Major depressive disorder | 1 | 34923209 | <https://doi.org/10.1016/j.euroneuro.2021.11.009> | NA |
| PainoldA.\_2019 | NA | NA | cross-sectional observational, not case-control | Bipolar disorder | 14 | 30051546 | <https://doi.org/10.1111/bdi.12682> | NA |

### Taxa table

``` r

createTaxonTable(my.dat.cond, n = 20) %>% 
  kbl() %>% 
  kable_styling()
```

| Taxon Name | Taxonomic Level | total_signatures | increased_signatures | decreased_signatures | Binomial Test pval |
|:---|:---|---:|---:|---:|---:|
| Streptococcus | genus | 9 | 2 | 7 | 0.18 |
| Coriobacteriaceae | family | 6 | 4 | 2 | 0.69 |
| Clostridiaceae | family | 6 | 2 | 4 | 0.69 |
| Vibrio | genus | 6 | 1 | 5 | 0.22 |
| Aerococcus | genus | 5 | 4 | 1 | 0.37 |
| Facklamia | genus | 5 | 1 | 4 | 0.37 |
| Enterococcus | genus | 5 | 1 | 4 | 0.37 |
| Faecalibacterium | genus | 5 | 1 | 4 | 0.37 |
| Psychrobacter | genus | 5 | 4 | 1 | 0.37 |
| Lactobacillaceae | family | 4 | 4 | 0 | 0.12 |
| Lactobacillus | genus | 4 | 4 | 0 | 0.12 |
| Streptococcaceae | family | 4 | 1 | 3 | 0.63 |
| Lachnospiraceae | family | 4 | 3 | 1 | 0.63 |
| Roseburia | genus | 4 | 0 | 4 | 0.12 |
| Peptostreptococcaceae | family | 4 | 1 | 3 | 0.63 |
| Megamonas | genus | 4 | 0 | 4 | 0.12 |
| Akkermansia | genus | 4 | 1 | 3 | 0.63 |
| Rothia | genus | 4 | 1 | 3 | 0.63 |
| Actinomycetota | phylum | 3 | 2 | 1 | 1.00 |
| Oscillospiraceae | family | 3 | 1 | 2 | 1.00 |

## Tables of only bipolar disorder

``` r

ind <- my.dat.cond[,"Condition"] %in% c("bipolar disorder")
my.dat.bd <- my.dat.cond[!ind,]
dim(my.dat.bd)
```

    ## [1] 50 51

``` r

table(my.dat.bd[,"Condition"])
```

    ## 
    ##          Bipolar disorder Major depressive disorder 
    ##                        18                        32

### Studies table

``` r

createStudyTable(my.dat.bd) %>% 
  kbl() %>% 
  kable_styling()
```

| Study Identifier | MaxCases | MaxControls | Study design | Condition | N_signatures | PMID | DOI | URL |
|:---|---:|---:|:---|:---|---:|:---|:---|:---|
| ChenJ.J.\_2020 | 45 | 44 | cross-sectional observational, not case-control | Major depressive disorder | 7 | 32040443 | <https://doi.org/10.18632/aging.102775> | <https://pubmed.ncbi.nlm.nih.gov/32040443/> |
| ChenY.\_2024 | 29 | 33 | case-control | Major depressive disorder | 4 | 38436041 | <https://doi.org/10.2147/NDT.S449224> | NA |
| CoelloK.\_2019 | 113 | 39 | cross-sectional observational, not case-control | Bipolar disorder | 1 | 30261302 | <https://doi.org/10.1016/j.bbi.2018.09.026> | NA |
| EvansS.J.\_2017 | 115 | 64 | case-control | Bipolar disorder | 1 | 27988330 | <https://doi.org/10.1016/j.jpsychires.2016.12.007> | NA |
| FlowersS.A.\_2017 | 49 | 68 | cross-sectional observational, not case-control | Bipolar disorder | 2 | 28035686 | <https://doi.org/10.1002/phar.1890> | NA |
| JiangY.\_2024 | 65 | 61 | prospective cohort | Major depressive disorder | 2 | 39695082 | <https://doi.org/10.1038/s41398-024-03208-z> | <https://www.nature.com/articles/s41398-024-03208-z> |
| LiY.\_2024 | NA | NA | laboratory experiment | Major depressive disorder | 13 | 38344423 | <https://doi.org/10.2147/NDT.S448940> | NA |
| LiuP.\_2024 | 85 | 85 | case-control,laboratory experiment | Major depressive disorder | 5 | 38589368 | <https://doi.org/10.1038/s41467-024-47273-w> | NA |
| Malan-MullerS.\_2022 | 20 | 2 | case-control | Major depressive disorder | 1 | 34923209 | <https://doi.org/10.1016/j.euroneuro.2021.11.009> | NA |
| PainoldA.\_2019 | NA | NA | cross-sectional observational, not case-control | Bipolar disorder | 14 | 30051546 | <https://doi.org/10.1111/bdi.12682> | NA |

### Taxa table

``` r

#createTaxonTable(my.dat.bd, n = 20) %>% 
#  kbl() %>% 
#   kable_styling()

createTaxonTable(my.dat.bd, n = 20) %>%  
  kbl() %>%
  kable_styling()
```

| Taxon Name | Taxonomic Level | total_signatures | increased_signatures | decreased_signatures | Binomial Test pval |
|:---|:---|---:|---:|---:|---:|
| Streptococcus | genus | 9 | 2 | 7 | 0.18 |
| Coriobacteriaceae | family | 6 | 4 | 2 | 0.69 |
| Clostridiaceae | family | 6 | 2 | 4 | 0.69 |
| Vibrio | genus | 6 | 1 | 5 | 0.22 |
| Aerococcus | genus | 5 | 4 | 1 | 0.37 |
| Facklamia | genus | 5 | 1 | 4 | 0.37 |
| Enterococcus | genus | 5 | 1 | 4 | 0.37 |
| Faecalibacterium | genus | 5 | 1 | 4 | 0.37 |
| Psychrobacter | genus | 5 | 4 | 1 | 0.37 |
| Lactobacillaceae | family | 4 | 4 | 0 | 0.12 |
| Lactobacillus | genus | 4 | 4 | 0 | 0.12 |
| Streptococcaceae | family | 4 | 1 | 3 | 0.63 |
| Lachnospiraceae | family | 4 | 3 | 1 | 0.63 |
| Roseburia | genus | 4 | 0 | 4 | 0.12 |
| Peptostreptococcaceae | family | 4 | 1 | 3 | 0.63 |
| Megamonas | genus | 4 | 0 | 4 | 0.12 |
| Akkermansia | genus | 4 | 1 | 3 | 0.63 |
| Rothia | genus | 4 | 1 | 3 | 0.63 |
| Actinomycetota | phylum | 3 | 2 | 1 | 1.00 |
| Oscillospiraceae | family | 3 | 1 | 2 | 1.00 |

## Tables only on major depressive disorder

``` r

ind <- my.dat.cond[,"Condition"] %in% c("major depressive disorder")
my.dat.mdd <- my.dat.cond[!ind,]
dim(my.dat.mdd)
```

    ## [1] 50 51

``` r

table(my.dat.mdd[,"Condition"])
```

    ## 
    ##          Bipolar disorder Major depressive disorder 
    ##                        18                        32

### Studies table

``` r

#createStudyTable(my.dat.mdd) %>% 
#  kbl() %>% 
 # kable_styling()

createStudyTable(my.dat.mdd) %>%
  kbl() %>%
  kable_styling()
```

| Study Identifier | MaxCases | MaxControls | Study design | Condition | N_signatures | PMID | DOI | URL |
|:---|---:|---:|:---|:---|---:|:---|:---|:---|
| ChenJ.J.\_2020 | 45 | 44 | cross-sectional observational, not case-control | Major depressive disorder | 7 | 32040443 | <https://doi.org/10.18632/aging.102775> | <https://pubmed.ncbi.nlm.nih.gov/32040443/> |
| ChenY.\_2024 | 29 | 33 | case-control | Major depressive disorder | 4 | 38436041 | <https://doi.org/10.2147/NDT.S449224> | NA |
| CoelloK.\_2019 | 113 | 39 | cross-sectional observational, not case-control | Bipolar disorder | 1 | 30261302 | <https://doi.org/10.1016/j.bbi.2018.09.026> | NA |
| EvansS.J.\_2017 | 115 | 64 | case-control | Bipolar disorder | 1 | 27988330 | <https://doi.org/10.1016/j.jpsychires.2016.12.007> | NA |
| FlowersS.A.\_2017 | 49 | 68 | cross-sectional observational, not case-control | Bipolar disorder | 2 | 28035686 | <https://doi.org/10.1002/phar.1890> | NA |
| JiangY.\_2024 | 65 | 61 | prospective cohort | Major depressive disorder | 2 | 39695082 | <https://doi.org/10.1038/s41398-024-03208-z> | <https://www.nature.com/articles/s41398-024-03208-z> |
| LiY.\_2024 | NA | NA | laboratory experiment | Major depressive disorder | 13 | 38344423 | <https://doi.org/10.2147/NDT.S448940> | NA |
| LiuP.\_2024 | 85 | 85 | case-control,laboratory experiment | Major depressive disorder | 5 | 38589368 | <https://doi.org/10.1038/s41467-024-47273-w> | NA |
| Malan-MullerS.\_2022 | 20 | 2 | case-control | Major depressive disorder | 1 | 34923209 | <https://doi.org/10.1016/j.euroneuro.2021.11.009> | NA |
| PainoldA.\_2019 | NA | NA | cross-sectional observational, not case-control | Bipolar disorder | 14 | 30051546 | <https://doi.org/10.1111/bdi.12682> | NA |

### Taxa table

``` r

createTaxonTable(my.dat.mdd, n = 20) %>% 
  kbl() %>% 
  kable_styling()
```

| Taxon Name | Taxonomic Level | total_signatures | increased_signatures | decreased_signatures | Binomial Test pval |
|:---|:---|---:|---:|---:|---:|
| Streptococcus | genus | 9 | 2 | 7 | 0.18 |
| Coriobacteriaceae | family | 6 | 4 | 2 | 0.69 |
| Clostridiaceae | family | 6 | 2 | 4 | 0.69 |
| Vibrio | genus | 6 | 1 | 5 | 0.22 |
| Aerococcus | genus | 5 | 4 | 1 | 0.37 |
| Facklamia | genus | 5 | 1 | 4 | 0.37 |
| Enterococcus | genus | 5 | 1 | 4 | 0.37 |
| Faecalibacterium | genus | 5 | 1 | 4 | 0.37 |
| Psychrobacter | genus | 5 | 4 | 1 | 0.37 |
| Lactobacillaceae | family | 4 | 4 | 0 | 0.12 |
| Lactobacillus | genus | 4 | 4 | 0 | 0.12 |
| Streptococcaceae | family | 4 | 1 | 3 | 0.63 |
| Lachnospiraceae | family | 4 | 3 | 1 | 0.63 |
| Roseburia | genus | 4 | 0 | 4 | 0.12 |
| Peptostreptococcaceae | family | 4 | 1 | 3 | 0.63 |
| Megamonas | genus | 4 | 0 | 4 | 0.12 |
| Akkermansia | genus | 4 | 1 | 3 | 0.63 |
| Rothia | genus | 4 | 1 | 3 | 0.63 |
| Actinomycetota | phylum | 3 | 2 | 1 | 1.00 |
| Oscillospiraceae | family | 3 | 1 | 2 | 1.00 |
