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
    ## ✔ ggplot2   4.0.4     ✔ tibble    3.3.1
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

    ## Using cached version from 2026-11-20 02:49:18

``` r

dim(dat)
```

    ## [1] 12086    57

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
| Gastric cancer                           |   13 |
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

    ## Using cached version from 2026-11-20 02:37:20

``` r

dat.bpd <- bugsigdbr::subsetByOntology(dat, column = "Condition", "bipolar disorder", efo)
dat.upd <- bugsigdbr::subsetByOntology(dat, column = "Condition", "major depressive disorder", efo)
my.dat.cond <- rbind(dat.bpd, dat.upd)
table(my.dat.cond[,"Condition"])
```

    ## 
    ##          Bipolar disorder Major depressive disorder       Unipolar depression 
    ##                        18                        40                        13

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
| JiangH.\_2015 | 29 | 30 | case-control | Unipolar depression | 4 | 25882912 | <https://doi.org/10.1016/j.bbi.2015.03.016> | NA |
| JiangY.\_2024 | 65 | 61 | prospective cohort | Major depressive disorder | 2 | 39695082 | <https://doi.org/10.1038/s41398-024-03208-z> | <https://www.nature.com/articles/s41398-024-03208-z> |
| LiY.\_2024 | NA | NA | laboratory experiment | Major depressive disorder | 13 | 38344423 | <https://doi.org/10.2147/NDT.S448940> | NA |
| LinP.\_2017 | 10 | 10 | case-control | Unipolar depression | 6 | 27741466 | <https://doi.org/10.1016/j.jad.2016.09.051> | NA |
| LingZ.\_2022 | 92 | 48 | case-control | Major depressive disorder | 4 | 36059521 | <https://doi.org/10.3389/fimmu.2022.964910> | <https://www.frontiersin.org/articles/10.3389/fimmu.2022.964910/full> |
| LiuF.\_2025 | 19 | 13 | cross-sectional observational, not case-control | Major depressive disorder | 2 | 41309831 | <https://doi.org/10.1038/s41598-025-26526-8> | <https://www.nature.com/articles/s41598-025-26526-8> |
| LiuP.\_2024 | 85 | 85 | case-control,laboratory experiment | Major depressive disorder | 5 | 38589368 | <https://doi.org/10.1038/s41467-024-47273-w> | NA |
| MadanA.\_2020 | NA | NA | prospective cohort | Unipolar depression | 3 | 32056780 | <https://doi.org/10.1016/j.jad.2019.12.020> | NA |
| Malan-MullerS.\_2022 | 20 | 2 | case-control | Major depressive disorder | 1 | 34923209 | <https://doi.org/10.1016/j.euroneuro.2021.11.009> | NA |
| PainoldA.\_2019 | NA | NA | cross-sectional observational, not case-control | Bipolar disorder | 14 | 30051546 | <https://doi.org/10.1111/bdi.12682> | NA |
| WingfieldB.\_2021 | 40 | 43 | case-control | Major depressive disorder | 2 | 34294835 | <https://doi.org/10.1038/s41598-021-94498-6> | <https://www.nature.com/articles/s41598-021-94498-6> |

### Taxa table

``` r

createTaxonTable(my.dat.cond, n = 20) %>% 
  kbl() %>% 
  kable_styling()
```

| Taxon Name | Taxonomic Level | total_signatures | increased_signatures | decreased_signatures | Binomial Test pval |
|:---|:---|---:|---:|---:|---:|
| Streptococcus | genus | 13 | 6 | 7 | 1.000 |
| Coriobacteriaceae | family | 8 | 6 | 2 | 0.290 |
| Faecalibacterium | genus | 8 | 1 | 7 | 0.070 |
| Peptostreptococcaceae | family | 8 | 4 | 4 | 1.000 |
| Clostridiaceae | family | 7 | 2 | 5 | 0.450 |
| Actinomycetota | phylum | 6 | 4 | 2 | 0.690 |
| Streptococcaceae | family | 6 | 3 | 3 | 1.000 |
| Lachnospiraceae | family | 6 | 5 | 1 | 0.220 |
| Bacteroidota | phylum | 6 | 1 | 5 | 0.220 |
| Prevotella | genus | 6 | 4 | 2 | 0.690 |
| Vibrio | genus | 6 | 1 | 5 | 0.220 |
| Bacillota | phylum | 5 | 5 | 0 | 0.062 |
| Aerococcus | genus | 5 | 4 | 1 | 0.370 |
| Facklamia | genus | 5 | 1 | 4 | 0.370 |
| Enterococcus | genus | 5 | 1 | 4 | 0.370 |
| Prevotellaceae | family | 5 | 3 | 2 | 1.000 |
| Psychrobacter | genus | 5 | 4 | 1 | 0.370 |
| Akkermansia | genus | 5 | 1 | 4 | 0.370 |
| Lactobacillaceae | family | 4 | 4 | 0 | 0.120 |
| Lactobacillus | genus | 4 | 4 | 0 | 0.120 |

## Tables of only bipolar disorder

``` r

ind <- my.dat.cond[,"Condition"] %in% c("bipolar disorder")
my.dat.bd <- my.dat.cond[!ind,]
dim(my.dat.bd)
```

    ## [1] 71 57

``` r

table(my.dat.bd[,"Condition"])
```

    ## 
    ##          Bipolar disorder Major depressive disorder       Unipolar depression 
    ##                        18                        40                        13

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
| JiangH.\_2015 | 29 | 30 | case-control | Unipolar depression | 4 | 25882912 | <https://doi.org/10.1016/j.bbi.2015.03.016> | NA |
| JiangY.\_2024 | 65 | 61 | prospective cohort | Major depressive disorder | 2 | 39695082 | <https://doi.org/10.1038/s41398-024-03208-z> | <https://www.nature.com/articles/s41398-024-03208-z> |
| LiY.\_2024 | NA | NA | laboratory experiment | Major depressive disorder | 13 | 38344423 | <https://doi.org/10.2147/NDT.S448940> | NA |
| LinP.\_2017 | 10 | 10 | case-control | Unipolar depression | 6 | 27741466 | <https://doi.org/10.1016/j.jad.2016.09.051> | NA |
| LingZ.\_2022 | 92 | 48 | case-control | Major depressive disorder | 4 | 36059521 | <https://doi.org/10.3389/fimmu.2022.964910> | <https://www.frontiersin.org/articles/10.3389/fimmu.2022.964910/full> |
| LiuF.\_2025 | 19 | 13 | cross-sectional observational, not case-control | Major depressive disorder | 2 | 41309831 | <https://doi.org/10.1038/s41598-025-26526-8> | <https://www.nature.com/articles/s41598-025-26526-8> |
| LiuP.\_2024 | 85 | 85 | case-control,laboratory experiment | Major depressive disorder | 5 | 38589368 | <https://doi.org/10.1038/s41467-024-47273-w> | NA |
| MadanA.\_2020 | NA | NA | prospective cohort | Unipolar depression | 3 | 32056780 | <https://doi.org/10.1016/j.jad.2019.12.020> | NA |
| Malan-MullerS.\_2022 | 20 | 2 | case-control | Major depressive disorder | 1 | 34923209 | <https://doi.org/10.1016/j.euroneuro.2021.11.009> | NA |
| PainoldA.\_2019 | NA | NA | cross-sectional observational, not case-control | Bipolar disorder | 14 | 30051546 | <https://doi.org/10.1111/bdi.12682> | NA |
| WingfieldB.\_2021 | 40 | 43 | case-control | Major depressive disorder | 2 | 34294835 | <https://doi.org/10.1038/s41598-021-94498-6> | <https://www.nature.com/articles/s41598-021-94498-6> |

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
| Streptococcus | genus | 13 | 6 | 7 | 1.000 |
| Coriobacteriaceae | family | 8 | 6 | 2 | 0.290 |
| Faecalibacterium | genus | 8 | 1 | 7 | 0.070 |
| Peptostreptococcaceae | family | 8 | 4 | 4 | 1.000 |
| Clostridiaceae | family | 7 | 2 | 5 | 0.450 |
| Actinomycetota | phylum | 6 | 4 | 2 | 0.690 |
| Streptococcaceae | family | 6 | 3 | 3 | 1.000 |
| Lachnospiraceae | family | 6 | 5 | 1 | 0.220 |
| Bacteroidota | phylum | 6 | 1 | 5 | 0.220 |
| Prevotella | genus | 6 | 4 | 2 | 0.690 |
| Vibrio | genus | 6 | 1 | 5 | 0.220 |
| Bacillota | phylum | 5 | 5 | 0 | 0.062 |
| Aerococcus | genus | 5 | 4 | 1 | 0.370 |
| Facklamia | genus | 5 | 1 | 4 | 0.370 |
| Enterococcus | genus | 5 | 1 | 4 | 0.370 |
| Prevotellaceae | family | 5 | 3 | 2 | 1.000 |
| Psychrobacter | genus | 5 | 4 | 1 | 0.370 |
| Akkermansia | genus | 5 | 1 | 4 | 0.370 |
| Lactobacillaceae | family | 4 | 4 | 0 | 0.120 |
| Lactobacillus | genus | 4 | 4 | 0 | 0.120 |

## Tables only on major depressive disorder

``` r

ind <- my.dat.cond[,"Condition"] %in% c("major depressive disorder")
my.dat.mdd <- my.dat.cond[!ind,]
dim(my.dat.mdd)
```

    ## [1] 71 57

``` r

table(my.dat.mdd[,"Condition"])
```

    ## 
    ##          Bipolar disorder Major depressive disorder       Unipolar depression 
    ##                        18                        40                        13

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
| JiangH.\_2015 | 29 | 30 | case-control | Unipolar depression | 4 | 25882912 | <https://doi.org/10.1016/j.bbi.2015.03.016> | NA |
| JiangY.\_2024 | 65 | 61 | prospective cohort | Major depressive disorder | 2 | 39695082 | <https://doi.org/10.1038/s41398-024-03208-z> | <https://www.nature.com/articles/s41398-024-03208-z> |
| LiY.\_2024 | NA | NA | laboratory experiment | Major depressive disorder | 13 | 38344423 | <https://doi.org/10.2147/NDT.S448940> | NA |
| LinP.\_2017 | 10 | 10 | case-control | Unipolar depression | 6 | 27741466 | <https://doi.org/10.1016/j.jad.2016.09.051> | NA |
| LingZ.\_2022 | 92 | 48 | case-control | Major depressive disorder | 4 | 36059521 | <https://doi.org/10.3389/fimmu.2022.964910> | <https://www.frontiersin.org/articles/10.3389/fimmu.2022.964910/full> |
| LiuF.\_2025 | 19 | 13 | cross-sectional observational, not case-control | Major depressive disorder | 2 | 41309831 | <https://doi.org/10.1038/s41598-025-26526-8> | <https://www.nature.com/articles/s41598-025-26526-8> |
| LiuP.\_2024 | 85 | 85 | case-control,laboratory experiment | Major depressive disorder | 5 | 38589368 | <https://doi.org/10.1038/s41467-024-47273-w> | NA |
| MadanA.\_2020 | NA | NA | prospective cohort | Unipolar depression | 3 | 32056780 | <https://doi.org/10.1016/j.jad.2019.12.020> | NA |
| Malan-MullerS.\_2022 | 20 | 2 | case-control | Major depressive disorder | 1 | 34923209 | <https://doi.org/10.1016/j.euroneuro.2021.11.009> | NA |
| PainoldA.\_2019 | NA | NA | cross-sectional observational, not case-control | Bipolar disorder | 14 | 30051546 | <https://doi.org/10.1111/bdi.12682> | NA |
| WingfieldB.\_2021 | 40 | 43 | case-control | Major depressive disorder | 2 | 34294835 | <https://doi.org/10.1038/s41598-021-94498-6> | <https://www.nature.com/articles/s41598-021-94498-6> |

### Taxa table

``` r

createTaxonTable(my.dat.mdd, n = 20) %>% 
  kbl() %>% 
  kable_styling()
```

| Taxon Name | Taxonomic Level | total_signatures | increased_signatures | decreased_signatures | Binomial Test pval |
|:---|:---|---:|---:|---:|---:|
| Streptococcus | genus | 13 | 6 | 7 | 1.000 |
| Coriobacteriaceae | family | 8 | 6 | 2 | 0.290 |
| Faecalibacterium | genus | 8 | 1 | 7 | 0.070 |
| Peptostreptococcaceae | family | 8 | 4 | 4 | 1.000 |
| Clostridiaceae | family | 7 | 2 | 5 | 0.450 |
| Actinomycetota | phylum | 6 | 4 | 2 | 0.690 |
| Streptococcaceae | family | 6 | 3 | 3 | 1.000 |
| Lachnospiraceae | family | 6 | 5 | 1 | 0.220 |
| Bacteroidota | phylum | 6 | 1 | 5 | 0.220 |
| Prevotella | genus | 6 | 4 | 2 | 0.690 |
| Vibrio | genus | 6 | 1 | 5 | 0.220 |
| Bacillota | phylum | 5 | 5 | 0 | 0.062 |
| Aerococcus | genus | 5 | 4 | 1 | 0.370 |
| Facklamia | genus | 5 | 1 | 4 | 0.370 |
| Enterococcus | genus | 5 | 1 | 4 | 0.370 |
| Prevotellaceae | family | 5 | 3 | 2 | 1.000 |
| Psychrobacter | genus | 5 | 4 | 1 | 0.370 |
| Akkermansia | genus | 5 | 1 | 4 | 0.370 |
| Lactobacillaceae | family | 4 | 4 | 0 | 0.120 |
| Lactobacillus | genus | 4 | 4 | 0 | 0.120 |
