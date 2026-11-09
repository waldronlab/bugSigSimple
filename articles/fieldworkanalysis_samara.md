# The Endometriosis-associated Microbiome

## Analysis objectives

1.  Import, recode, and subset data from bugsigdb.org  
2.  Create a table of studies
3.  Calculate the frequency of appearance of each taxa in independent
    signatures and identify the most frequently reported taxa
4.  Estimate the probability of the most frequently identified taxa
    occuring by chance

## Making sure packages are installed

``` r

install.packages(c("devtools", "tidyverse", "kableExtra"))
devtools::install_github("waldronlab/BugSigDBStats")
devtools::install_github("waldronlab/bugsigdbr")
```

## Load and subset data

``` r

library(bugSigSimple)
library(BugSigDBStats)
library(bugsigdbr)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.2.1     ✔ readr     2.2.0
    ## ✔ forcats   1.0.1     ✔ stringr   1.6.0
    ## ✔ ggplot2   4.0.3     ✔ tibble    3.3.1
    ## ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
    ## ✔ purrr     1.2.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r

library(stringr)
library(kableExtra)
```

    ## 
    ## Attaching package: 'kableExtra'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     group_rows

``` r

library(dplyr)
```

``` r

library(bugSigSimple)
dat <- bugsigdbr::importBugSigDB(cache= TRUE)
```

    ## Using cached version from 2026-11-09 16:32:08

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

Subsetting only studies done on humans - 6 of 17 studies were excluded
because they were either done on mice or not statistically reliable
(Cregger et. al)

``` r

included.pmid <-
  c(
    "30778155",
    "32192080",
    "31087436",
    "26901400",
    "33839907",
    "32046455",
    "33925708",
    "32299442",
    "33313185",
    "34268384",
    "33660232"
  )
subset.dat <- dat[dat$PMID %in% included.pmid, ]
```

## All studies

Summary of studies and most frequent taxa increased and decreased in
endometriosis patients for all studies

``` r

bugSigSimple::createStudyTable(subset.dat) %>%
  kbl() %>%
  kable_paper("hover", full_width = FALSE) %>%
  scroll_box(width = "100%", height = "400px")
```

| Study Identifier | MaxCases | MaxControls | Study design | Condition | N_signatures | PMID | DOI | URL |
|:---|---:|---:|:---|:---|---:|:---|:---|:---|
| AkiyamaK.\_2019 | 30 | 39 | cross-sectional observational, not case-control | Endometriosis | 1 | 31087436 | <https://doi.org/10.1111/aji.13147> | NA |
| AtaB.\_2019 | 14 | 14 | prospective cohort | Endometriosis | 2 | 30778155 | <https://doi.org/10.1038/s41598-019-39700-6> | <https://pubmed.ncbi.nlm.nih.gov/30778155/> |
| ChaoX.\_2021 | 37 | 66 | case-control | Endometriosis | 4 | 34268384 | <https://doi.org/10.21037/atm-20-4586> | NA |
| ChenS.\_2020 | 14 | 120 | case-control | Adenomyosis,Endometriosis | 3 | 33313185 | <https://doi.org/10.21037/atm-20-1309> | NA |
| KhanK.N.\_2016 | 32 | 32 | case-control | Endometriosis | 5 | 26901400 | <https://doi.org/10.1016/j.ejogrb.2016.01.040> | NA |
| LeeS.R.\_2021 | 45 | 45 | case-control | Endometriosis | 2 | 33925708 | <https://doi.org/10.3390/ijms22094608> | NA |
| PerrottaA.R.\_2020 | 8 | 12 | cross-sectional observational, not case-control | Endometriosis | 2 | 32046455 | <https://doi.org/10.1007/s43032-019-00113-5> | NA |
| ShanJ.\_2021 | 12 | 12 | case-control | Endometriosis | 2 | 33839907 | <https://doi.org/10.1007/s00404-021-06057-z> | NA |
| SvenssonA.\_2021 | 66 | 198 | case-control | Endometriosis | 4 | 33660232 | <https://doi.org/10.1007/s43032-021-00506-5> | NA |
| WeiW.\_2020 | 36 | 14 | cross-sectional observational, not case-control | Endometriosis | 2 | 32299442 | <https://doi.org/10.1186/s12941-020-00356-0> | NA |

``` r

BugSigDBStats::getMostFrequentTaxa(subset.dat, n=30) %>%
  kbl() %>%
  kable_paper("hover", full_width = FALSE) %>%
  scroll_box(height = "300px")
```

| Var1    | Freq |
|:--------|-----:|
| 1301    |    6 |
| 1578    |    4 |
| 543     |    4 |
| 84107   |    4 |
| 119852  |    3 |
| 1300    |    3 |
| 135621  |    3 |
| 171549  |    3 |
| 1940338 |    3 |
| 28050   |    3 |
| 39948   |    3 |
| 469     |    3 |
| 838     |    3 |
| 906     |    3 |
| 1350    |    2 |
| 1380    |    2 |
| 1508657 |    2 |
| 1582879 |    2 |
| 1654    |    2 |
| 1663    |    2 |
| 168808  |    2 |
| 1743    |    2 |
| 194     |    2 |
| 216851  |    2 |
| 2701    |    2 |
| 286     |    2 |
| 33038   |    2 |
| 33958   |    2 |
| 419014  |    2 |
| 72294   |    2 |

``` r

BugSigDBStats::getMostFrequentTaxa(subset.dat, direction="UP") %>%
  kbl() %>%
  kable_paper("hover", full_width = FALSE) %>%
  scroll_box(height = "300px")
```

| Var1    | Freq |
|:--------|-----:|
| 1301    |    5 |
| 543     |    4 |
| 119852  |    3 |
| 1300    |    3 |
| 135621  |    3 |
| 1578    |    3 |
| 1940338 |    3 |
| 469     |    3 |
| 84107   |    3 |
| 1350    |    2 |

``` r

BugSigDBStats::getMostFrequentTaxa(subset.dat, direction="DOWN") %>%
  kbl() %>%
  kable_paper("hover", full_width = FALSE) %>%
  scroll_box(height = "300px")
```

| Var1    | Freq |
|:--------|-----:|
| 168808  |    2 |
| 171549  |    2 |
| 28050   |    2 |
| 33958   |    2 |
| 906     |    2 |
| 1185407 |    1 |
| 1301    |    1 |
| 1452070 |    1 |
| 1508657 |    1 |
| 1578    |    1 |

## Excluding feces samples

Summary of studies and most frequent taxa in only samples from female
reproductive tract, excluding feces samples

``` r

included.pmid2 <-
  c(
    "30778155",
    "32192080",
    "31087436",
    "26901400",
    "32046455",
    "33925708",
    "32299442",
    "33313185",
    "34268384"
  )
subset.dat2 <- dat[dat$PMID %in% included.pmid2, ]

reproductive_sigs <-
  subset.dat2[which(subset.dat2$`Body site` != "feces" |
                      is.na(subset.dat2$`Body site`)), ]

bugSigSimple::createStudyTable(reproductive_sigs) %>%
  kbl() %>%
  kable_paper("hover", full_width = FALSE) %>%
  scroll_box(width = "100%", height = "400px")
```

| Study Identifier | MaxCases | MaxControls | Study design | Condition | N_signatures | PMID | DOI | URL |
|:---|---:|---:|:---|:---|---:|:---|:---|:---|
| AkiyamaK.\_2019 | 30 | 39 | cross-sectional observational, not case-control | Endometriosis | 1 | 31087436 | <https://doi.org/10.1111/aji.13147> | NA |
| AtaB.\_2019 | 14 | 14 | prospective cohort | Endometriosis | 2 | 30778155 | <https://doi.org/10.1038/s41598-019-39700-6> | <https://pubmed.ncbi.nlm.nih.gov/30778155/> |
| ChaoX.\_2021 | 37 | 66 | case-control | Endometriosis | 4 | 34268384 | <https://doi.org/10.21037/atm-20-4586> | NA |
| ChenS.\_2020 | 14 | 120 | case-control | Adenomyosis,Endometriosis | 3 | 33313185 | <https://doi.org/10.21037/atm-20-1309> | NA |
| KhanK.N.\_2016 | 32 | 32 | case-control | Endometriosis | 5 | 26901400 | <https://doi.org/10.1016/j.ejogrb.2016.01.040> | NA |
| LeeS.R.\_2021 | 45 | 45 | case-control | Endometriosis | 2 | 33925708 | <https://doi.org/10.3390/ijms22094608> | NA |
| PerrottaA.R.\_2020 | 8 | 12 | cross-sectional observational, not case-control | Endometriosis | 2 | 32046455 | <https://doi.org/10.1007/s43032-019-00113-5> | NA |
| WeiW.\_2020 | 36 | 14 | cross-sectional observational, not case-control | Endometriosis | 2 | 32299442 | <https://doi.org/10.1186/s12941-020-00356-0> | NA |

``` r

allfreqs <- bugSigSimple::createTaxonTable(reproductive_sigs, n = 20) %>% #could change number
  arrange(I(decreased_signatures - increased_signatures))
incfreqs <- filter(allfreqs, I(increased_signatures - decreased_signatures) > 0)
decfreqs <- filter(allfreqs, I(increased_signatures - decreased_signatures) < 0)
kableExtra::kbl(allfreqs) %>%
  kable_paper("hover", full_width = FALSE)
```

| Taxon Name | Taxonomic Level | total_signatures | increased_signatures | decreased_signatures | Binomial Test pval |
|:---|:---|---:|---:|---:|---:|
| Streptococcus | genus | 5 | 4 | 1 | 0.37 |
| Coriobacteriaceae | family | 3 | 3 | 0 | 0.25 |
| Streptococcaceae | family | 3 | 3 | 0 | 0.25 |
| Enterobacteriaceae | family | 3 | 3 | 0 | 0.25 |
| Escherichia/Shigella sp. | species | 3 | 3 | 0 | 0.25 |
| Acinetobacter | genus | 3 | 3 | 0 | 0.25 |
| Pseudomonadaceae | family | 3 | 3 | 0 | 0.25 |
| Alloscardovia | genus | 2 | 2 | 0 | 0.50 |
| Arthrobacter | genus | 2 | 2 | 0 | 0.50 |
| Atopobium | genus | 2 | 2 | 0 | 0.50 |
| Staphylococcaceae | family | 2 | 2 | 0 | 0.50 |
| Faecalibacterium | genus | 2 | 2 | 0 | 0.50 |
| Ezakiella | genus | 2 | 2 | 0 | 0.50 |
| Lactobacillus | genus | 3 | 2 | 1 | 1.00 |
| Dialister | genus | 3 | 2 | 1 | 1.00 |
| Prevotella | genus | 3 | 2 | 1 | 1.00 |
| Gardnerella | genus | 2 | 1 | 1 | 1.00 |
| Propionibacterium | genus | 2 | 1 | 1 | 1.00 |
| Megasphaera | genus | 3 | 1 | 2 | 1.00 |
| Lactobacillaceae | family | 2 | 0 | 2 | 0.50 |

### Load HMP data from curatedMetagenomicData

``` r

suppressPackageStartupMessages(library(curatedMetagenomicData))
# hmp_vagina_metadata <- filter(sampleMetadata, study_name == "HMP_2012" & body_site == "vagina")
all_healthy_adult_vagina_metadata <- filter(sampleMetadata, body_site == "vagina" & disease == "healthy" & age_category == "adult")
se <-
  curatedMetagenomicData::returnSamples(all_healthy_adult_vagina_metadata, dataType = "relative_abundance", counts = FALSE)
```

Get matrices of species and genus relative abundance

``` r

allranks <- mia::splitByRanks(se)
#species_relab <- t(assay(allranks[["species"]]))
genus_relab <- t(assay(allranks[["genus"]]))
#family_relab <- t(assay(allranks[["family"]]))
```

Centered log-ratio transformation on genus level relative abundance:

``` r

genus_relab_clr <- Hotelling::clr(genus_relab + 1)
```

``` r

df <- dplyr::select(data.frame(genus_relab_clr), c("Streptococcus", "Lactobacillus"))
df$age <- se$age
df$BMI <- se$BMI
fit <- lm(Streptococcus ~ Lactobacillus + age + BMI, data = df)  
summary(fit)
```

    ## 
    ## Call:
    ## lm(formula = Streptococcus ~ Lactobacillus + age + BMI, data = df)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.6728 -0.1875 -0.0264  0.0893  3.6729 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   -0.386413   0.436694  -0.885 0.378987    
    ## Lactobacillus -0.078095   0.047933  -1.629 0.107344    
    ## age            0.035738   0.010073   3.548 0.000665 ***
    ## BMI           -0.009094   0.014295  -0.636 0.526581    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4857 on 77 degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## Multiple R-squared:  0.1799, Adjusted R-squared:  0.1479 
    ## F-statistic: 5.629 on 3 and 77 DF,  p-value: 0.00153

``` r

summary(lm(Streptococcus ~ Lactobacillus, data = df))
```

    ## 
    ## Call:
    ## lm(formula = Streptococcus ~ Lactobacillus, data = df)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.6006 -0.0609 -0.0462 -0.0440  4.1909 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)    0.55246    0.20193   2.736  0.00745 **
    ## Lactobacillus -0.11409    0.04751  -2.401  0.01832 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5642 on 93 degrees of freedom
    ## Multiple R-squared:  0.05839,    Adjusted R-squared:  0.04826 
    ## F-statistic: 5.767 on 1 and 93 DF,  p-value: 0.01832

``` r

wilcox.test(df$Streptococcus, df$Lactobacillus)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  df$Streptococcus and df$Lactobacillus
    ## W = 213.5, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0

``` r

t.test(df$Streptococcus, df$Lactobacillus)
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  df$Streptococcus and df$Lactobacillus
    ## t = -28.666, df = 133.93, p-value < 2.2e-16
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -4.258800 -3.709054
    ## sample estimates:
    ## mean of x mean of y 
    ## 0.0879008 4.0718280

``` r

secount <-
  curatedMetagenomicData::returnSamples(all_healthy_adult_vagina_metadata, dataType = "relative_abundance", counts = TRUE)
```

    ## 
    ## $`2021-03-31.FerrettiP_2018.relative_abundance`
    ## dropping rows without rowTree matches:
    ##   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Coriobacteriaceae|g__Collinsella|s__Collinsella_stercoris
    ##   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Coriobacteriaceae|g__Enorma|s__[Collinsella]_massiliensis
    ##   k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Bacillales_unclassified|g__Gemella|s__Gemella_bergeri
    ##   k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Carnobacteriaceae|g__Granulicatella|s__Granulicatella_elegans
    ##   k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcus|s__Ruminococcus_champanellensis
    ##   k__Bacteria|p__Firmicutes|c__Erysipelotrichia|o__Erysipelotrichales|f__Erysipelotrichaceae|g__Bulleidia|s__Bulleidia_extructa
    ##   k__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Sutterellaceae|g__Sutterella|s__Sutterella_parvirubra
    ##   k__Bacteria|p__Spirochaetes|c__Spirochaetia|o__Leptospirales|f__Leptospiraceae|g__Leptospira|s__Leptospira_hartskeerlii
    ##   k__Bacteria|p__Synergistetes|c__Synergistia|o__Synergistales|f__Synergistaceae|g__Cloacibacillus|s__Cloacibacillus_evryensis

    ## $`2021-03-31.HMP_2012.relative_abundance`
    ## dropping rows without rowTree matches:
    ##   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Atopobiaceae|g__Olsenella|s__Olsenella_profusa
    ##   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Coriobacteriaceae|g__Collinsella|s__Collinsella_stercoris
    ##   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Coriobacteriaceae|g__Enorma|s__[Collinsella]_massiliensis
    ##   k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Bacillales_unclassified|g__Gemella|s__Gemella_bergeri
    ##   k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Carnobacteriaceae|g__Granulicatella|s__Granulicatella_elegans
    ##   k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcus|s__Ruminococcus_champanellensis
    ##   k__Bacteria|p__Firmicutes|c__Erysipelotrichia|o__Erysipelotrichales|f__Erysipelotrichaceae|g__Bulleidia|s__Bulleidia_extructa
    ##   k__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Sutterellaceae|g__Sutterella|s__Sutterella_parvirubra
    ##   k__Bacteria|p__Synergistetes|c__Synergistia|o__Synergistales|f__Synergistaceae|g__Cloacibacillus|s__Cloacibacillus_evryensis

    ## $`2021-03-31.WampachL_2018.relative_abundance`
    ## dropping rows without rowTree matches:
    ##   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Coriobacteriaceae|g__Collinsella|s__Collinsella_stercoris
    ##   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Coriobacteriaceae|g__Enorma|s__[Collinsella]_massiliensis
    ##   k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Carnobacteriaceae|g__Granulicatella|s__Granulicatella_elegans
    ##   k__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Sutterellaceae|g__Sutterella|s__Sutterella_parvirubra
    ##   k__Bacteria|p__Synergistetes|c__Synergistia|o__Synergistales|f__Synergistaceae|g__Cloacibacillus|s__Cloacibacillus_evryensis

``` r

allranks <- mia::splitByRanks(secount)
genus_relab_count <- t(assay(allranks[["genus"]]))
df <- dplyr::select(data.frame(genus_relab_count), c("Streptococcus", "Lactobacillus"))
df$age <- secount$age
df$BMI <- secount$BMI
```

``` r

library(pscl)
```

    ## Classes and Methods for R originally developed in the
    ## Political Science Computational Laboratory
    ## Department of Political Science
    ## Stanford University (2002-2015),
    ## by and under the direction of Simon Jackman.
    ## hurdle and zeroinfl functions by Achim Zeileis.

``` r

fit.ZInegbin <- pscl::zeroinfl(Streptococcus ~ Lactobacillus + age + BMI | 1,
                data = df,
                dist = "negbin")
```

## Create a data frame for regression analysis from healthy controls for log-linear model

``` r

library(curatedMetagenomicData)
library(dplyr)
library(DT)
library(purrr)

taxon_strep <- c("1884" = "Streptococcus")
taxon_lact <- c("1578" = "Lactobacillus")

target_sample <- "vagina"
sample_metadata <- sampleMetadata %>% 
    filter(
        grepl(target_sample, body_site) | grepl(target_sample, body_subsite)
    ) %>% 
    discard(~all(is.na(.x))) %>% 
    as_tibble()
dim(sample_metadata)
```

    ## [1] 96 26

``` r

unique(sample_metadata$study_name)
```

    ## [1] "FerrettiP_2018" "HMP_2012"       "WampachL_2018"

``` r

tse <- returnSamples(
    sampleMetadata = sample_metadata,
    dataType = "relative_abundance",
    counts = FALSE,
    rownames = "short"
)
```

    ## 
    ## $`2021-03-31.FerrettiP_2018.relative_abundance`
    ## dropping rows without rowTree matches:
    ##   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Coriobacteriaceae|g__Collinsella|s__Collinsella_stercoris
    ##   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Coriobacteriaceae|g__Enorma|s__[Collinsella]_massiliensis
    ##   k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Bacillales_unclassified|g__Gemella|s__Gemella_bergeri
    ##   k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Carnobacteriaceae|g__Granulicatella|s__Granulicatella_elegans
    ##   k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcus|s__Ruminococcus_champanellensis
    ##   k__Bacteria|p__Firmicutes|c__Erysipelotrichia|o__Erysipelotrichales|f__Erysipelotrichaceae|g__Bulleidia|s__Bulleidia_extructa
    ##   k__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Sutterellaceae|g__Sutterella|s__Sutterella_parvirubra
    ##   k__Bacteria|p__Spirochaetes|c__Spirochaetia|o__Leptospirales|f__Leptospiraceae|g__Leptospira|s__Leptospira_hartskeerlii
    ##   k__Bacteria|p__Synergistetes|c__Synergistia|o__Synergistales|f__Synergistaceae|g__Cloacibacillus|s__Cloacibacillus_evryensis

    ## Duplicated labels were made unique.

    ## $`2021-03-31.HMP_2012.relative_abundance`
    ## dropping rows without rowTree matches:
    ##   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Atopobiaceae|g__Olsenella|s__Olsenella_profusa
    ##   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Coriobacteriaceae|g__Collinsella|s__Collinsella_stercoris
    ##   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Coriobacteriaceae|g__Enorma|s__[Collinsella]_massiliensis
    ##   k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Bacillales_unclassified|g__Gemella|s__Gemella_bergeri
    ##   k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Carnobacteriaceae|g__Granulicatella|s__Granulicatella_elegans
    ##   k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcus|s__Ruminococcus_champanellensis
    ##   k__Bacteria|p__Firmicutes|c__Erysipelotrichia|o__Erysipelotrichales|f__Erysipelotrichaceae|g__Bulleidia|s__Bulleidia_extructa
    ##   k__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Sutterellaceae|g__Sutterella|s__Sutterella_parvirubra
    ##   k__Bacteria|p__Synergistetes|c__Synergistia|o__Synergistales|f__Synergistaceae|g__Cloacibacillus|s__Cloacibacillus_evryensis

    ## Duplicated labels were made unique.

    ## $`2021-03-31.WampachL_2018.relative_abundance`
    ## dropping rows without rowTree matches:
    ##   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Coriobacteriaceae|g__Collinsella|s__Collinsella_stercoris
    ##   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Coriobacteriaceae|g__Enorma|s__[Collinsella]_massiliensis
    ##   k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Carnobacteriaceae|g__Granulicatella|s__Granulicatella_elegans
    ##   k__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Sutterellaceae|g__Sutterella|s__Sutterella_parvirubra
    ##   k__Bacteria|p__Synergistetes|c__Synergistia|o__Synergistales|f__Synergistaceae|g__Cloacibacillus|s__Cloacibacillus_evryensis

    ## Duplicated labels were made unique.

    ## Warning in convertNode(tree = tree, node = labN, use.alias = FALSE, message =
    ## FALSE): Multiple nodes are found to have the same label.

``` r

tse
```

    ## class: TreeSummarizedExperiment 
    ## dim: 908 96 
    ## metadata(0):
    ## assays(1): relative_abundance
    ## rownames(908): species:Escherichia coli species:Klebsiella
    ##   michiganensis ... species:Pediococcus pentosaceus species:Rahnella
    ##   aquatilis
    ## rowData names(7): superkingdom phylum ... genus species
    ## colnames(96): CA_C10001MS2003VA_t0M15 CA_C10003MS2031VA_t0M15 ...
    ##   MV_C118 MV_C120
    ## colData names(25): study_name subject_id ... family_role body_subsite
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(0):
    ## rowLinks: a LinkDataFrame (908 rows)
    ## rowTree: 1 phylo tree(s) (10430 leaves)
    ## colLinks: NULL
    ## colTree: NULL

### Comparing mean to variance of the outcome variable

``` r

library(ggplot2)
mean(df$Streptococcus)
```

    ## [1] 4049.558

``` r

var(df$Streptococcus)
```

    ## [1] 509068831

``` r

fraczeroes <- nrow(df[df$Streptococcus == 0,])/nrow(df)
fraczeroes
```

    ## [1] 0.8736842

\#Fit log-linear models

``` r

fit.pois <- glm(Streptococcus ~ Lactobacillus + age + BMI,
                data = df,
                family = poisson(link = "log"))

summary (fit.pois)
```

    ## 
    ## Call:
    ## glm(formula = Streptococcus ~ Lactobacillus + age + BMI, family = poisson(link = "log"), 
    ##     data = df)
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   -5.313e-01  1.873e-02  -28.36   <2e-16 ***
    ## Lactobacillus  5.141e-08  9.828e-10   52.32   <2e-16 ***
    ## age            3.599e-01  4.136e-04  870.25   <2e-16 ***
    ## BMI           -1.573e-01  6.453e-04 -243.77   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 1929049  on 80  degrees of freedom
    ## Residual deviance:  632598  on 77  degrees of freedom
    ##   (14 observations deleted due to missingness)
    ## AIC: 632704
    ## 
    ## Number of Fisher Scoring iterations: 20

``` r

library(pscl)
fit.ZIpois <- zeroinfl(Streptococcus ~ Lactobacillus + age + BMI | 1,
                data = df,
                dist = "poisson")

summary (fit.ZIpois)
```

    ## 
    ## Call:
    ## zeroinfl(formula = Streptococcus ~ Lactobacillus + age + BMI | 1, data = df, 
    ##     dist = "poisson")
    ## 
    ## Pearson residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.3753 -0.3753 -0.3752 -0.3750 27.8168 
    ## 
    ## Count model coefficients (poisson with log link):
    ##                 Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)    6.169e+00        NaN     NaN      NaN
    ## Lactobacillus  1.321e-07        NaN     NaN      NaN
    ## age            2.420e-01        NaN     NaN      NaN
    ## BMI           -2.116e-01        NaN     NaN      NaN
    ## 
    ## Zero-inflation model coefficients (binomial with logit link):
    ##             Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)     1.96        NaN     NaN      NaN
    ## 
    ## Number of iterations in BFGS optimization: 20 
    ## Log-likelihood: -1.02e+05 on 5 Df

``` r

library(MASS)
```

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r

fit.negbin <-
  glm.nb(
    Streptococcus ~ Lactobacillus + age + BMI,
    data = df,
    control = glm.control(
      epsilon = 1e-8,
      maxit = 50,
      trace = FALSE
    )
  )
```

    ## Warning: glm.fit: algorithm did not converge

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning: glm.fit: algorithm did not converge

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning: glm.fit: algorithm did not converge

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning: glm.fit: algorithm did not converge

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning: glm.fit: algorithm did not converge

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning: glm.fit: algorithm did not converge

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning: glm.fit: algorithm did not converge

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning: glm.fit: algorithm did not converge

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning: glm.fit: algorithm did not converge

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning: glm.fit: algorithm did not converge

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning: glm.fit: algorithm did not converge

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning: glm.fit: algorithm did not converge

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning: glm.fit: algorithm did not converge

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in glm.nb(Streptococcus ~ Lactobacillus + age + BMI, data = df, :
    ## alternation limit reached

``` r

library(broom)
tidy(fit.negbin)
```

    ## # A tibble: 4 × 5
    ##   term               estimate   std.error statistic     p.value
    ##   <chr>                 <dbl>       <dbl>     <dbl>       <dbl>
    ## 1 (Intercept)   -34.4         8.11            -4.24 0.0000224  
    ## 2 Lactobacillus   0.000000889 0.000000196      4.54 0.00000570 
    ## 3 age             1.09        0.207            5.30 0.000000119
    ## 4 BMI             0.297       0.277            1.07 0.282

``` r

fit.ZInegbin <- zeroinfl(Streptococcus ~ Lactobacillus + age + BMI | 1,
                data = df,
                dist = "negbin")

summary(fit.ZInegbin)
```

    ## 
    ## Call:
    ## zeroinfl(formula = Streptococcus ~ Lactobacillus + age + BMI | 1, data = df, 
    ##     dist = "negbin")
    ## 
    ## Pearson residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.1923 -0.1923 -0.1923 -0.1923  7.1121 
    ## 
    ## Count model coefficients (negbin with log link):
    ##                 Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)    1.917e+01        NaN     NaN      NaN
    ## Lactobacillus -2.725e-07        NaN     NaN      NaN
    ## age            3.274e-02        NaN     NaN      NaN
    ## BMI           -4.033e-01        NaN     NaN      NaN
    ## Log(theta)    -9.618e-01        NaN     NaN      NaN
    ## 
    ## Zero-inflation model coefficients (binomial with logit link):
    ##             Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)     1.91        NaN     NaN      NaN
    ## 
    ## Theta = 0.3822 
    ## Number of iterations in BFGS optimization: 27 
    ## Log-likelihood: -131.9 on 6 Df
