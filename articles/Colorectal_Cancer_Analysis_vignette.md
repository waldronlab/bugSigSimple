# Exploratory Analysis of Colorectal Cancer (CRC) Microbial Signatures

## Introduction

This vignette demonstrates an exploratory analysis of microbial
signatures associated with colorectal cancer (CRC) using data obtained
from BugSigDB and functions provided by the bugSigSimple package.

## Analysis objectives

1.  Downloading the data from BugSigDB
2.  Filtering of the data to select relevant signatures and remove
    similar/duplicate signatures
3.  Restricting to fecal samples (the body site with the most
    consistent, human, non-tissue evidence)
4.  Description of the studies
5.  Description of the taxa
6.  Basic statistical tests
7.  Monte Carlo simulation to calculate critical N
8.  Heatmaps and dendrograms

## Making sure packages are installed

``` r

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("devtools", "tidyverse", "kableExtra"))
BiocManager::install(c("waldronlab/bugSigSimple", "waldronlab/BugSigDBStats", "waldronlab/bugsigdbr"))
```

## Download BugSigDB data

``` r

full.dat <- bugsigdbr::importBugSigDB(cache = FALSE)
```

``` r

dim(full.dat)
```

    ## [1] 12086    57

The BugSigDB database was downloaded using the
[`importBugSigDB()`](https://rdrr.io/pkg/bugsigdbr/man/importBugSigDB.html)
function from the **bugsigdbr** package. The result contains all curated
microbial signatures together with study metadata such as disease, body
site, sequencing platform, and abundance direction.

## Remove empty signatures

``` r

ind1 <- lengths(full.dat[["MetaPhlAn taxon names"]]) > 0
ind2 <- lengths(full.dat[["NCBI Taxonomy IDs"]]) > 0

dat <- full.dat[ind1 & ind2, ]

dim(dat)
```

    ## [1] 12086    57

Signatures without annotated taxa were removed to ensure that only
complete microbial signatures were included in downstream analyses.

## Select CRC-related signatures

``` r

pattern <- paste(
  "Colorectal cancer",
  "Colorectal Neoplasms",
  "Colorectal carcinoma",
  "Colorectal adenocarcinoma",
  "Colonic Neoplasms",
  "Rectal Neoplasms",
  "Colorectal adenoma",
  "Colorectal polyp",
  "colon cancer",
  "rectal cancer",
  "CRC",
  sep = "|"
)

crc.dat <- dat %>%
  filter(
    grepl(pattern, Condition, ignore.case = TRUE) |
      grepl(pattern, `Group 0 name`, ignore.case = TRUE) |
      grepl(pattern, `Group 1 name`, ignore.case = TRUE)
  ) %>%
  filter(
    grepl(pattern, Condition, ignore.case = TRUE)
  )

nrow(crc.dat)
```

    ## [1] 487

``` r

crc.dat %>%
  count(Condition, sort = TRUE)
```

    ##                                               Condition   n
    ## 1                                     Colorectal cancer 374
    ## 2                                    Colorectal adenoma  54
    ## 3                                  Colorectal carcinoma  13
    ## 4                          Metastatic colorectal cancer  12
    ## 5                             Colorectal adenocarcinoma  11
    ## 6                   Colorectal cancer,Diabetes mellitus   5
    ## 7                       Nutraceutical,Colorectal cancer   4
    ## 8                       Colorectal cancer,Nutraceutical   3
    ## 9                  Colorectal adenoma,Colorectal cancer   2
    ## 10 Colorectal adenoma,Non-alcoholic fatty liver disease   2
    ## 11                          Colorectal cancer cell line   2
    ## 12            Colorectal cancer,Peripheral nerve injury   2
    ## 13 Non-alcoholic fatty liver disease,Colorectal adenoma   2
    ## 14                Colorectal cancer,Response to placebo   1

To capture the complete body of evidence related to colorectal cancer
(CRC), signatures were identified using both the primary disease
terminology (colorectal cancer / neoplasms) and related terms covering
adenomas, polyps, and site-specific variants (colon/rectal cancer).
Disease-related terms were searched across the Condition, Group 0 name,
and Group 1 name fields to recover studies curated under different
naming conventions. Finally, only records whose primary condition
corresponded to a CRC-related disease were retained.

## Remove similar signatures

``` r

crc.filtered <- crc.dat %>%
  distinct()

nrow(crc.filtered)
```

    ## [1] 487

Exact duplicate records were removed using
[`distinct()`](https://dplyr.tidyverse.org/reference/distinct.html) to
reduce redundancy before downstream analyses. This step retains unique
CRC-associated microbial signatures while preserving the overall study
representation.

## Restrict to fecal samples

``` r

crc.filtered %>%
  count(`Body site`, sort = TRUE) %>%
  knitr::kable()
```

| Body site                                       |   n |
|:------------------------------------------------|----:|
| Feces                                           | 311 |
| Intestine                                       |  44 |
| Colon                                           |  23 |
| Colorectal mucosa                               |  23 |
| Colorectum                                      |  16 |
| Saliva                                          |  12 |
| Oral cavity                                     |  10 |
| feces                                           |  10 |
| Ascending colon,Colorectal mucosa,Sigmoid colon |   6 |
| Sigmoid colon                                   |   5 |
| Ascending colon                                 |   4 |
| Colon,Rectum,Sigmoid colon                      |   4 |
| Colonic mucosa                                  |   4 |
| Intestinal mucosa                               |   4 |
| Ascending colon,Colorectal mucosa               |   2 |
| Ascending colon,Sigmoid colon,Colorectal mucosa |   2 |
| Colorectal mucosa,Sigmoid colon                 |   2 |
| Dorsum of tongue                                |   2 |
| Mucosa of ascending colon                       |   2 |
| Rectum                                          |   1 |

``` r

# Feces is the dominant, human, non-tissue body site in the CRC signatures.
# One study (PMID 26397137) used a mouse model rather than human subjects and
# is excluded so that the signature set reflects human fecal microbiome findings.
crc.filtered <- crc.filtered %>%
  filter(`Body site` == "Feces") %>%
  filter(PMID != 26397137)

nrow(crc.filtered)
```

    ## [1] 309

CRC signatures in BugSigDB were curated from several body sites,
including feces, mucosal tissue, and oral samples, as well as from a
small number of mouse studies. To keep the downstream comparisons
consistent, the analysis was restricted to fecal signatures from human
studies, which represent the largest and most directly comparable body
of evidence for the gut microbiome in CRC.

## Study characteristics

``` r

study.summary <-
  crc.filtered %>%
  count(`Body site`, `Sequencing type`)

knitr::kable(study.summary)
```

| Body site | Sequencing type |   n |
|:----------|:----------------|----:|
| Feces     | 16S             | 209 |
| Feces     | WMS             | 100 |

The selected CRC fecal signatures span multiple studies and sequencing
approaches. Summarizing study characteristics provides context for
interpreting downstream analyses and highlights the diversity of
experimental designs represented in BugSigDB.

## Taxonomic summary

``` r

taxa.table <- createTaxonTable(crc.filtered)

head(taxa.table)
```

    ## # A tibble: 6 × 6
    ##   `Taxon Name`       `Taxonomic Level` total_signatures increased_signatures
    ##   <chr>              <chr>                        <int>                <int>
    ## 1 Parvimonas         genus                           29                   29
    ## 2 Peptostreptococcus genus                           27                   27
    ## 3 Bacteroides        genus                           27                   10
    ## 4 Clostridium        genus                           26                    8
    ## 5 Lachnospiraceae    family                          26                    6
    ## 6 Ruminococcus       genus                           24                   11
    ## # ℹ 2 more variables: decreased_signatures <int>, `Binomial Test pval` <dbl>

``` r

knitr::kable(head(taxa.table, 10))
```

| Taxon Name | Taxonomic Level | total_signatures | increased_signatures | decreased_signatures | Binomial Test pval |
|:---|:---|---:|---:|---:|---:|
| Parvimonas | genus | 29 | 29 | 0 | 0.0e+00 |
| Peptostreptococcus | genus | 27 | 27 | 0 | 0.0e+00 |
| Bacteroides | genus | 27 | 10 | 17 | 2.5e-01 |
| Clostridium | genus | 26 | 8 | 18 | 7.6e-02 |
| Lachnospiraceae | family | 26 | 6 | 20 | 9.4e-03 |
| Ruminococcus | genus | 24 | 11 | 13 | 8.4e-01 |
| Fusobacterium | genus | 23 | 19 | 4 | 2.6e-03 |
| Fusobacterium nucleatum | species | 23 | 22 | 1 | 5.7e-06 |
| Bacillota | phylum | 22 | 13 | 9 | 5.2e-01 |
| Prevotella | genus | 22 | 11 | 11 | 1.0e+00 |

The taxonomic summary identifies microorganisms that recur across
CRC-associated signatures. Frequently reported taxa may represent
microbial features that are consistently observed across independent CRC
studies.

## Split by sequencing type

``` r

crc.wms <- crc.filtered %>% filter(`Sequencing type` == "WMS")
crc.16s <- crc.filtered %>% filter(`Sequencing type` == "16S")

cat("WMS  - signatures:", nrow(crc.wms), "| studies:", n_distinct(crc.wms$PMID), "\n")
```

    ## WMS  - signatures: 100 | studies: 18

``` r

cat("16S  - signatures:", nrow(crc.16s), "| studies:", n_distinct(crc.16s$PMID), "\n")
```

    ## 16S  - signatures: 209 | studies: 57

Splitting by sequencing platform makes it possible to check later
whether whole-metagenome shotgun (WMS) and 16S rRNA studies agree on the
direction of the taxa they report.

## Long-format taxon table

``` r

sigs.feces <- bugsigdbr::getSignatures(crc.filtered, tax.id.type = "taxname")

taxon.long <- purrr::map_dfr(seq_along(sigs.feces), function(i) {
  taxa <- sigs.feces[[i]]
  if (length(taxa) == 0) return(NULL)
  tibble::tibble(sig_id = i, taxon = taxa)
}) %>%
  left_join(
    crc.filtered %>%
      mutate(sig_id = row_number()) %>%
      select(sig_id, PMID, Condition, `Sequencing type`, `Abundance in Group 1`),
    by = "sig_id"
  ) %>%
  rename(direction = `Abundance in Group 1`)

nrow(taxon.long)
```

    ## [1] 2986

Each signature was expanded into one row per taxon and joined back to
its study metadata (PMID, condition, sequencing type, and reported
direction). This long-format table underlies all of the frequency,
concordance, enrichment, and co-occurrence analyses below.

## Taxon frequency across studies

``` r

taxon.freq <- taxon.long %>%
  distinct(PMID, taxon, direction) %>%
  count(taxon, direction, name = "n_studies") %>%
  group_by(taxon) %>%
  mutate(total_studies = sum(n_studies)) %>%
  ungroup() %>%
  arrange(desc(total_studies))

n_distinct(taxon.freq$taxon)
```

    ## [1] 907

``` r

taxon.freq %>%
  distinct(taxon, total_studies) %>%
  slice_max(total_studies, n = 10) %>%
  knitr::kable()
```

| taxon            | total_studies |
|:-----------------|--------------:|
| Bacteroides      |            20 |
| Lachnospiraceae  |            20 |
| Fusobacterium    |            18 |
| Ruminococcus     |            18 |
| Clostridium      |            17 |
| Faecalibacterium |            17 |
| Bacillota        |            16 |
| Coprococcus      |            16 |
| Prevotella       |            16 |
| Roseburia        |            16 |

The most frequently reported taxa give a first impression of which
microbes are consistently associated with CRC across independent
studies, before accounting for direction of change or statistical
significance.

## Core consensus signature

``` r

concordance <- taxon.long %>%
  distinct(PMID, taxon, direction) %>%
  group_by(taxon) %>%
  summarise(
    n_studies       = n_distinct(PMID),
    n_increased     = sum(direction == "increased"),
    n_decreased     = sum(direction == "decreased"),
    concordance_pct = round(max(n_increased, n_decreased) / n_studies * 100, 1),
    majority_dir    = if_else(n_increased >= n_decreased, "increased", "decreased"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_studies), desc(concordance_pct))

# Formal threshold for a "core" signature: reported in at least 3 independent
# studies, with at least 75% of those studies agreeing on direction
core.sig <- concordance %>%
  filter(n_studies >= 3, concordance_pct >= 75) %>%
  arrange(majority_dir, desc(n_studies))

nrow(core.sig)
```

    ## [1] 135

``` r

knitr::kable(core.sig)
```

| taxon | n_studies | n_increased | n_decreased | concordance_pct | majority_dir |
|:---|---:|---:|---:|---:|:---|
| Lachnospiraceae | 18 | 6 | 14 | 77.8 | decreased |
| Faecalibacterium | 16 | 3 | 14 | 87.5 | decreased |
| Clostridium | 15 | 5 | 12 | 80.0 | decreased |
| Lachnospira | 15 | 3 | 12 | 80.0 | decreased |
| Clostridia | 12 | 3 | 9 | 75.0 | decreased |
| Bifidobacterium | 11 | 3 | 9 | 81.8 | decreased |
| Faecalibacterium prausnitzii | 9 | 1 | 9 | 100.0 | decreased |
| Blautia | 9 | 1 | 8 | 88.9 | decreased |
| Anaerobutyricum hallii | 8 | 1 | 7 | 87.5 | decreased |
| Lachnospira eligens | 8 | 1 | 7 | 87.5 | decreased |
| Subdoligranulum | 8 | 1 | 7 | 87.5 | decreased |
| Actinomycetota | 7 | 1 | 6 | 85.7 | decreased |
| Romboutsia | 7 | 1 | 6 | 85.7 | decreased |
| Roseburia intestinalis | 6 | 0 | 6 | 100.0 | decreased |
| Bifidobacteriaceae | 6 | 1 | 5 | 83.3 | decreased |
| Bifidobacteriales | 6 | 1 | 5 | 83.3 | decreased |
| Segatella copri | 6 | 2 | 5 | 83.3 | decreased |
| Bifidobacterium animalis | 5 | 0 | 5 | 100.0 | decreased |
| Bifidobacterium longum | 5 | 0 | 5 | 100.0 | decreased |
| Staphylococcus | 5 | 1 | 4 | 80.0 | decreased |
| Agathobacter | 4 | 0 | 4 | 100.0 | decreased |
| Agathobacter rectalis | 4 | 0 | 4 | 100.0 | decreased |
| Anaerostipes hadrus | 4 | 1 | 4 | 100.0 | decreased |
| Bifidobacterium adolescentis | 4 | 0 | 4 | 100.0 | decreased |
| Cyanobacteriota | 4 | 1 | 4 | 100.0 | decreased |
| Fusicatenibacter | 4 | 0 | 4 | 100.0 | decreased |
| Pseudocoprococcus catus | 4 | 0 | 4 | 100.0 | decreased |
| Tyzzerella | 4 | 0 | 4 | 100.0 | decreased |
| Anaerotruncus | 4 | 1 | 3 | 75.0 | decreased |
| Bacteroides intestinalis | 4 | 1 | 3 | 75.0 | decreased |
| Bacteroides thetaiotaomicron | 4 | 1 | 3 | 75.0 | decreased |
| Bifidobacterium catenulatum | 4 | 1 | 3 | 75.0 | decreased |
| Clostridiales bacterium | 4 | 2 | 3 | 75.0 | decreased |
| Eubacterium coprostanoligenes | 4 | 1 | 3 | 75.0 | decreased |
| Leptotrichia | 4 | 1 | 3 | 75.0 | decreased |
| Micrococcaceae | 4 | 1 | 3 | 75.0 | decreased |
| Neisseria | 4 | 1 | 3 | 75.0 | decreased |
| Allocoprococcus comes | 3 | 0 | 3 | 100.0 | decreased |
| Blautia obeum | 3 | 0 | 3 | 100.0 | decreased |
| Coprobacillus | 3 | 0 | 3 | 100.0 | decreased |
| Roseburia faecis | 3 | 1 | 3 | 100.0 | decreased |
| Streptococcus salivarius | 3 | 0 | 3 | 100.0 | decreased |
| Streptococcus thermophilus | 3 | 0 | 3 | 100.0 | decreased |
| Fusobacterium | 18 | 14 | 4 | 77.8 | increased |
| Parvimonas | 15 | 15 | 0 | 100.0 | increased |
| Peptostreptococcus | 15 | 15 | 0 | 100.0 | increased |
| Porphyromonas | 15 | 12 | 3 | 80.0 | increased |
| Fusobacterium nucleatum | 14 | 13 | 1 | 92.9 | increased |
| Enterobacteriaceae | 14 | 11 | 3 | 78.6 | increased |
| Pseudomonadota | 12 | 9 | 3 | 75.0 | increased |
| Parvimonas micra | 11 | 11 | 0 | 100.0 | increased |
| Fusobacteriota | 11 | 9 | 2 | 81.8 | increased |
| Parabacteroides | 11 | 9 | 3 | 81.8 | increased |
| Gemella morbillorum | 10 | 10 | 0 | 100.0 | increased |
| Peptostreptococcus stomatis | 10 | 10 | 0 | 100.0 | increased |
| Solobacterium moorei | 9 | 9 | 0 | 100.0 | increased |
| Bacteroides fragilis | 9 | 8 | 1 | 88.9 | increased |
| Eggerthella | 9 | 8 | 1 | 88.9 | increased |
| Lactobacillus | 9 | 8 | 2 | 88.9 | increased |
| \[Clostridium\] symbiosum | 9 | 8 | 2 | 88.9 | increased |
| Enterococcus | 9 | 7 | 3 | 77.8 | increased |
| Fusobacteriia | 9 | 7 | 2 | 77.8 | increased |
| Gemella | 8 | 8 | 0 | 100.0 | increased |
| Akkermansia | 8 | 7 | 1 | 87.5 | increased |
| Escherichia | 8 | 7 | 1 | 87.5 | increased |
| Akkermansia muciniphila | 8 | 6 | 2 | 75.0 | increased |
| Enterobacterales | 8 | 6 | 2 | 75.0 | increased |
| Escherichia/Shigella sp. | 8 | 6 | 3 | 75.0 | increased |
| Fusobacteriaceae | 8 | 6 | 2 | 75.0 | increased |
| Fusobacteriales | 8 | 6 | 2 | 75.0 | increased |
| Lactobacillales | 8 | 6 | 3 | 75.0 | increased |
| Escherichia coli | 7 | 7 | 0 | 100.0 | increased |
| Peptostreptococcus anaerobius | 7 | 7 | 0 | 100.0 | increased |
| Bacilli | 7 | 6 | 1 | 85.7 | increased |
| Porphyromonas asaccharolytica | 7 | 6 | 1 | 85.7 | increased |
| Lactobacillaceae | 6 | 6 | 1 | 100.0 | increased |
| Prevotella intermedia | 6 | 6 | 0 | 100.0 | increased |
| Verrucomicrobiota | 6 | 6 | 1 | 100.0 | increased |
| Gammaproteobacteria | 6 | 5 | 1 | 83.3 | increased |
| Bilophila wadsworthia | 5 | 5 | 0 | 100.0 | increased |
| Flavonifractor plautii | 5 | 5 | 0 | 100.0 | increased |
| Hungatella hathewayi | 5 | 5 | 1 | 100.0 | increased |
| Odoribacter splanchnicus | 5 | 5 | 0 | 100.0 | increased |
| Streptococcus anginosus | 5 | 5 | 1 | 100.0 | increased |
| Dialister pneumosintes | 5 | 4 | 1 | 80.0 | increased |
| Morganella | 5 | 4 | 1 | 80.0 | increased |
| Oscillospira | 5 | 4 | 1 | 80.0 | increased |
| Paraprevotella | 5 | 4 | 1 | 80.0 | increased |
| Proteus | 5 | 4 | 1 | 80.0 | increased |
| Selenomonadales | 5 | 4 | 1 | 80.0 | increased |
| Shigella | 5 | 4 | 1 | 80.0 | increased |
| Alistipes onderdonkii | 4 | 4 | 0 | 100.0 | increased |
| Bacteroidaceae | 4 | 4 | 0 | 100.0 | increased |
| Bacteroides salyersiae | 4 | 4 | 0 | 100.0 | increased |
| Hungatella | 4 | 4 | 0 | 100.0 | increased |
| Porphyromonas uenonis | 4 | 4 | 0 | 100.0 | increased |
| Verrucomicrobiales | 4 | 4 | 0 | 100.0 | increased |
| Verrucomicrobiia | 4 | 4 | 0 | 100.0 | increased |
| uncultured Oscillospiraceae bacterium | 4 | 4 | 0 | 100.0 | increased |
| Acinetobacter | 4 | 3 | 1 | 75.0 | increased |
| Anaerostipes caccae | 4 | 3 | 1 | 75.0 | increased |
| Atopobium | 4 | 3 | 1 | 75.0 | increased |
| Bacteroides ovatus | 4 | 3 | 1 | 75.0 | increased |
| Christensenella | 4 | 3 | 1 | 75.0 | increased |
| Corynebacteriaceae | 4 | 3 | 1 | 75.0 | increased |
| Desulfovibrionaceae | 4 | 3 | 1 | 75.0 | increased |
| Eggerthella lenta | 4 | 3 | 1 | 75.0 | increased |
| Enterobacter | 4 | 3 | 1 | 75.0 | increased |
| Negativicutes | 4 | 3 | 1 | 75.0 | increased |
| Phocaeicola vulgatus | 4 | 3 | 1 | 75.0 | increased |
| Ralstonia | 4 | 3 | 1 | 75.0 | increased |
| Streptococcus lutetiensis | 4 | 3 | 2 | 75.0 | increased |
| Thomasclavelia | 4 | 3 | 1 | 75.0 | increased |
| Abiotrophia | 3 | 3 | 0 | 100.0 | increased |
| Actinomycetaceae | 3 | 3 | 1 | 100.0 | increased |
| Aeromonas | 3 | 3 | 0 | 100.0 | increased |
| Alistipes finegoldii | 3 | 3 | 0 | 100.0 | increased |
| Alloprevotella tannerae | 3 | 3 | 0 | 100.0 | increased |
| Anaerococcus vaginalis | 3 | 3 | 0 | 100.0 | increased |
| Bacteroides eggerthii | 3 | 3 | 0 | 100.0 | increased |
| Bifidobacterium breve | 3 | 3 | 1 | 100.0 | increased |
| Comamonadaceae | 3 | 3 | 0 | 100.0 | increased |
| Desulfovibrio desulfuricans | 3 | 3 | 0 | 100.0 | increased |
| Eisenbergiella | 3 | 3 | 0 | 100.0 | increased |
| Klebsiella oxytoca | 3 | 3 | 0 | 100.0 | increased |
| Lachnospiraceae bacterium 3_1_57FAA_CT1 | 3 | 3 | 0 | 100.0 | increased |
| Methanobrevibacter smithii | 3 | 3 | 0 | 100.0 | increased |
| Olsenella | 3 | 3 | 0 | 100.0 | increased |
| Parabacteroides distasonis | 3 | 3 | 1 | 100.0 | increased |
| Parabacteroides merdae | 3 | 3 | 0 | 100.0 | increased |
| Selenomonas sputigena | 3 | 3 | 0 | 100.0 | increased |
| Solobacterium | 3 | 3 | 0 | 100.0 | increased |
| Streptococcus oralis | 3 | 3 | 0 | 100.0 | increased |
| Tissierellia | 3 | 3 | 0 | 100.0 | increased |
| Verrucomicrobiaceae | 3 | 3 | 0 | 100.0 | increased |

``` r

core.sig %>%
  mutate(
    direction_label = if_else(majority_dir == "increased",
                               "Increased in CRC", "Decreased in CRC"),
    # Wrap long taxon names onto a second line instead of letting them
    # overlap or run off the plot area
    taxon_label = paste0(str_wrap(taxon, width = 28), " (", n_studies, ")")
  ) %>%
  ggplot(aes(x = reorder(taxon_label, n_studies), y = n_studies, fill = direction_label)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(concordance_pct, "%")), hjust = -0.15, size = 3.3) +
  coord_flip(clip = "off") +
  facet_wrap(~direction_label, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("Increased in CRC" = "#B22222",
                                "Decreased in CRC" = "#006400")) +
  # Extra headroom on the right so the "%" labels never get clipped
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(
    title = "Core CRC Consensus Signature",
    subtitle = "Taxa in >=3 studies with >=75% direction concordance",
    x = NULL, y = "Number of independent studies", fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10, lineheight = 0.9),
    strip.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
  )
```

![](Colorectal_Cancer_Analysis_vignette_files/figure-html/core-signature-plot-1.png)

Requiring at least three independent studies and at least 75%
directional agreement is a common, interpretable threshold for calling a
taxon “core”: it rules out one-off findings while still tolerating some
disagreement between studies that differ in cohort, region, or platform.

## WMS vs 16S concordance

``` r

dir.by.platform <- taxon.long %>%
  filter(`Sequencing type` %in% c("WMS", "16S")) %>%
  distinct(PMID, taxon, direction, `Sequencing type`) %>%
  count(taxon, `Sequencing type`, direction, name = "n_studies") %>%
  group_by(taxon, `Sequencing type`) %>%
  slice_max(n_studies, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(taxon, `Sequencing type`, majority_dir = direction)

platform.compare <- dir.by.platform %>%
  pivot_wider(names_from = `Sequencing type`, values_from = majority_dir) %>%
  filter(!is.na(WMS), !is.na(`16S`)) %>%
  rename(dir_WMS = WMS, dir_16S = `16S`) %>%
  mutate(agree = dir_WMS == dir_16S)

nrow(platform.compare)
```

    ## [1] 222

``` r

sum(platform.compare$agree)
```

    ## [1] 142

``` r

round(100 * mean(platform.compare$agree), 1)
```

    ## [1] 64

``` r

knitr::kable(platform.compare)
```

| taxon                                     | dir_16S   | dir_WMS   | agree |
|:------------------------------------------|:----------|:----------|:------|
| Acidaminococcaceae                        | decreased | increased | FALSE |
| Acidaminococcus                           | decreased | decreased | TRUE  |
| Acinetobacter                             | increased | decreased | FALSE |
| Actinomyces                               | increased | decreased | FALSE |
| Actinomycetaceae                          | increased | increased | TRUE  |
| Actinomycetota                            | decreased | decreased | TRUE  |
| Aggregatibacter                           | increased | decreased | FALSE |
| Akkermansia muciniphila                   | increased | decreased | FALSE |
| Alistipes                                 | increased | decreased | FALSE |
| Alistipes finegoldii                      | increased | increased | TRUE  |
| Alistipes indistinctus                    | increased | increased | TRUE  |
| Alistipes onderdonkii                     | increased | increased | TRUE  |
| Alistipes putredinis                      | decreased | decreased | TRUE  |
| Alistipes senegalensis                    | decreased | increased | FALSE |
| Allocoprococcus comes                     | decreased | decreased | TRUE  |
| Alloprevotella tannerae                   | increased | increased | TRUE  |
| Anaerobutyricum hallii                    | decreased | decreased | TRUE  |
| Anaerococcus vaginalis                    | increased | increased | TRUE  |
| Anaerostipes                              | increased | decreased | FALSE |
| Anaerostipes caccae                       | increased | decreased | FALSE |
| Anaerostipes hadrus                       | decreased | decreased | TRUE  |
| Anaerotruncus                             | decreased | decreased | TRUE  |
| Anaerotruncus sp.                         | decreased | increased | FALSE |
| Atopobium                                 | increased | decreased | FALSE |
| Bacilli                                   | increased | increased | TRUE  |
| Bacillota                                 | increased | decreased | FALSE |
| Bacillus                                  | increased | decreased | FALSE |
| Bacteroides                               | decreased | decreased | TRUE  |
| Bacteroides caccae                        | increased | decreased | FALSE |
| Bacteroides eggerthii                     | increased | increased | TRUE  |
| Bacteroides finegoldii                    | decreased | decreased | TRUE  |
| Bacteroides fragilis                      | increased | increased | TRUE  |
| Bacteroides intestinalis                  | decreased | decreased | TRUE  |
| Bacteroides nordii                        | decreased | increased | FALSE |
| Bacteroides ovatus                        | increased | increased | TRUE  |
| Bacteroides salyersiae                    | increased | increased | TRUE  |
| Bacteroides sp.                           | decreased | increased | FALSE |
| Bacteroides stercoris                     | decreased | increased | FALSE |
| Bacteroides thetaiotaomicron              | decreased | decreased | TRUE  |
| Bacteroides uniformis                     | decreased | decreased | TRUE  |
| Bacteroidota                              | increased | decreased | FALSE |
| Bifidobacterium                           | decreased | decreased | TRUE  |
| Bifidobacterium adolescentis              | decreased | decreased | TRUE  |
| Bifidobacterium animalis                  | decreased | decreased | TRUE  |
| Bifidobacterium breve                     | increased | increased | TRUE  |
| Bifidobacterium catenulatum               | increased | decreased | FALSE |
| Bifidobacterium dentium                   | decreased | decreased | TRUE  |
| Bifidobacterium longum                    | decreased | decreased | TRUE  |
| Bifidobacterium pseudocatenulatum         | decreased | decreased | TRUE  |
| Bilophila                                 | increased | decreased | FALSE |
| Bilophila wadsworthia                     | increased | increased | TRUE  |
| Blautia                                   | decreased | decreased | TRUE  |
| Blautia obeum                             | decreased | decreased | TRUE  |
| Blautia producta                          | decreased | increased | FALSE |
| Burkholderiales                           | increased | decreased | FALSE |
| Butyrivibrio                              | increased | decreased | FALSE |
| Campylobacter rectus                      | increased | increased | TRUE  |
| Caryophanales                             | decreased | increased | FALSE |
| Christensenella                           | increased | decreased | FALSE |
| Christensenellaceae                       | decreased | decreased | TRUE  |
| Chryseobacterium                          | increased | decreased | FALSE |
| Citrobacter                               | decreased | decreased | TRUE  |
| Citrobacter freundii                      | increased | increased | TRUE  |
| Clostridia                                | decreased | decreased | TRUE  |
| Clostridiaceae                            | decreased | decreased | TRUE  |
| Clostridiales Family XVII. Incertae Sedis | increased | decreased | FALSE |
| Clostridiales bacterium 1_7_47FAA         | increased | increased | TRUE  |
| Clostridium                               | decreased | decreased | TRUE  |
| Clostridium sp.                           | decreased | decreased | TRUE  |
| Collinsella                               | decreased | decreased | TRUE  |
| Collinsella aerofaciens                   | increased | decreased | FALSE |
| Colwelliaceae                             | increased | decreased | FALSE |
| Coprobacter                               | increased | decreased | FALSE |
| Coprococcus                               | decreased | decreased | TRUE  |
| Desulfovibrio                             | decreased | decreased | TRUE  |
| Desulfovibrio desulfuricans               | increased | increased | TRUE  |
| Dialister                                 | decreased | decreased | TRUE  |
| Dialister pneumosintes                    | increased | increased | TRUE  |
| Dorea                                     | decreased | decreased | TRUE  |
| Dorea formicigenerans                     | decreased | increased | FALSE |
| Eggerthella                               | increased | increased | TRUE  |
| Eggerthella lenta                         | increased | decreased | FALSE |
| Enterobacter                              | increased | decreased | FALSE |
| Enterobacterales                          | increased | increased | TRUE  |
| Enterobacteriaceae                        | increased | decreased | FALSE |
| Enterococcus                              | increased | decreased | FALSE |
| Enterorhabdus                             | decreased | decreased | TRUE  |
| Erysipelotrichaceae                       | increased | decreased | FALSE |
| Erysipelotrichales                        | decreased | decreased | TRUE  |
| Erysipelotrichia                          | increased | decreased | FALSE |
| Escherichia                               | increased | increased | TRUE  |
| Escherichia coli                          | increased | increased | TRUE  |
| Eubacteriaceae                            | increased | decreased | FALSE |
| Eubacteriales                             | decreased | decreased | TRUE  |
| Eubacterium                               | increased | decreased | FALSE |
| Eubacterium ventriosum                    | decreased | decreased | TRUE  |
| Facklamia                                 | increased | decreased | FALSE |
| Faecalibacterium                          | decreased | decreased | TRUE  |
| Faecalibacterium prausnitzii              | decreased | decreased | TRUE  |
| Flavonifractor plautii                    | increased | increased | TRUE  |
| Fusicatenibacter                          | decreased | decreased | TRUE  |
| Fusobacteriia                             | increased | increased | TRUE  |
| Fusobacteriota                            | increased | increased | TRUE  |
| Fusobacterium                             | increased | decreased | FALSE |
| Fusobacterium nucleatum                   | increased | increased | TRUE  |
| Gemella                                   | increased | increased | TRUE  |
| Gemella morbillorum                       | increased | increased | TRUE  |
| Gordonibacter                             | increased | decreased | FALSE |
| Granulicatella adiacens                   | increased | increased | TRUE  |
| Haemophilus                               | decreased | decreased | TRUE  |
| Holdemanella biformis                     | increased | decreased | FALSE |
| Holdemania                                | decreased | decreased | TRUE  |
| Hungatella                                | increased | increased | TRUE  |
| Hungatella hathewayi                      | increased | increased | TRUE  |
| Intestinimonas                            | decreased | decreased | TRUE  |
| Klebsiella                                | decreased | decreased | TRUE  |
| Klebsiella oxytoca                        | increased | increased | TRUE  |
| Lachnoanaerobaculum                       | increased | decreased | FALSE |
| Lachnoclostridium                         | increased | decreased | FALSE |
| Lachnospira                               | decreased | decreased | TRUE  |
| Lachnospira eligens                       | decreased | decreased | TRUE  |
| Lachnospira pectinoschiza                 | decreased | decreased | TRUE  |
| Lachnospiraceae                           | decreased | decreased | TRUE  |
| Lachnospiraceae bacterium                 | decreased | increased | FALSE |
| Lachnospiraceae bacterium 3_1_57FAA_CT1   | increased | increased | TRUE  |
| Lactobacillaceae                          | increased | decreased | FALSE |
| Lactobacillales                           | increased | decreased | FALSE |
| Lactobacillus                             | increased | decreased | FALSE |
| Lactococcus                               | increased | decreased | FALSE |
| Leptotrichia                              | decreased | increased | FALSE |
| Leyella stercorea                         | decreased | increased | FALSE |
| Mediterraneibacter torques                | increased | increased | TRUE  |
| Megasphaera                               | increased | decreased | FALSE |
| Megasphaera elsdenii                      | decreased | decreased | TRUE  |
| Methanobrevibacter smithii                | increased | increased | TRUE  |
| Mitsuokella                               | decreased | decreased | TRUE  |
| Negativibacillus                          | decreased | decreased | TRUE  |
| Negativicutes                             | increased | increased | TRUE  |
| Neisseria                                 | decreased | decreased | TRUE  |
| Odoribacter                               | decreased | increased | FALSE |
| Olsenella                                 | increased | increased | TRUE  |
| Oscillibacter                             | decreased | decreased | TRUE  |
| Oscillospiraceae                          | decreased | decreased | TRUE  |
| Parabacteroides                           | increased | increased | TRUE  |
| Parabacteroides distasonis                | decreased | increased | FALSE |
| Parabacteroides goldsteinii               | decreased | decreased | TRUE  |
| Parabacteroides merdae                    | increased | increased | TRUE  |
| Parahaemophilus parainfluenzae            | decreased | decreased | TRUE  |
| Paraprevotella                            | increased | decreased | FALSE |
| Parasutterella                            | decreased | decreased | TRUE  |
| Parvimonas                                | increased | increased | TRUE  |
| Parvimonas micra                          | increased | increased | TRUE  |
| Parvimonas sp.                            | increased | increased | TRUE  |
| Pasteurellaceae                           | increased | decreased | FALSE |
| Pasteurellales                            | decreased | decreased | TRUE  |
| Peptoniphilaceae                          | increased | increased | TRUE  |
| Peptostreptococcaceae                     | increased | decreased | FALSE |
| Peptostreptococcus                        | increased | increased | TRUE  |
| Peptostreptococcus anaerobius             | increased | increased | TRUE  |
| Peptostreptococcus sp.                    | increased | increased | TRUE  |
| Peptostreptococcus stomatis               | increased | increased | TRUE  |
| Phocaeicola dorei                         | decreased | increased | FALSE |
| Phocaeicola vulgatus                      | decreased | increased | FALSE |
| Porphyromonadaceae                        | increased | increased | TRUE  |
| Porphyromonas                             | increased | increased | TRUE  |
| Porphyromonas asaccharolytica             | increased | increased | TRUE  |
| Porphyromonas somerae                     | increased | increased | TRUE  |
| Porphyromonas uenonis                     | increased | increased | TRUE  |
| Prevotella                                | increased | decreased | FALSE |
| Prevotella intermedia                     | increased | increased | TRUE  |
| Pseudobutyrivibrio                        | decreased | decreased | TRUE  |
| Pseudocoprococcus catus                   | decreased | decreased | TRUE  |
| Pseudoflavonifractor                      | decreased | increased | FALSE |
| Pseudoflavonifractor capillosus           | decreased | increased | FALSE |
| Pseudomonadota                            | increased | increased | TRUE  |
| Pyramidobacter                            | increased | decreased | FALSE |
| Ralstonia                                 | increased | decreased | FALSE |
| Robinsoniella                             | decreased | decreased | TRUE  |
| Romboutsia                                | decreased | decreased | TRUE  |
| Roseburia                                 | decreased | decreased | TRUE  |
| Roseburia intestinalis                    | decreased | decreased | TRUE  |
| Roseburia inulinivorans                   | decreased | increased | FALSE |
| Ruminococcus                              | increased | decreased | FALSE |
| Ruminococcus bromii                       | decreased | decreased | TRUE  |
| Ruminococcus sp. 5_1_39BFAA               | decreased | decreased | TRUE  |
| Salmonella                                | increased | decreased | FALSE |
| Segatella copri                           | decreased | decreased | TRUE  |
| Segatella oris                            | decreased | increased | FALSE |
| Selenomonadales                           | increased | increased | TRUE  |
| Selenomonas                               | increased | decreased | FALSE |
| Shigella                                  | increased | decreased | FALSE |
| Slackia sp.                               | increased | increased | TRUE  |
| Solobacterium                             | increased | increased | TRUE  |
| Solobacterium moorei                      | increased | increased | TRUE  |
| Streptococcaceae                          | increased | decreased | FALSE |
| Streptococcus                             | increased | decreased | FALSE |
| Streptococcus anginosus                   | increased | increased | TRUE  |
| Streptococcus lutetiensis                 | increased | decreased | FALSE |
| Streptococcus mutans                      | decreased | decreased | TRUE  |
| Streptococcus oralis                      | increased | increased | TRUE  |
| Streptococcus parasanguinis               | increased | increased | TRUE  |
| Streptococcus salivarius                  | decreased | decreased | TRUE  |
| Streptococcus thermophilus                | decreased | decreased | TRUE  |
| Subdoligranulum                           | decreased | decreased | TRUE  |
| Subdoligranulum sp. 4_3_54A2FAA           | increased | increased | TRUE  |
| Sutterellaceae                            | increased | decreased | FALSE |
| Synergistes sp. 3_1_syn1                  | increased | increased | TRUE  |
| Thomasclavelia                            | increased | decreased | FALSE |
| Tissierellales                            | increased | increased | TRUE  |
| Tissierellia                              | increased | increased | TRUE  |
| Turicibacter                              | increased | decreased | FALSE |
| Veillonella                               | decreased | decreased | TRUE  |
| Veillonella atypica                       | decreased | increased | FALSE |
| Veillonella dispar                        | decreased | decreased | TRUE  |
| Veillonella parvula                       | decreased | decreased | TRUE  |
| Veillonellaceae                           | increased | increased | TRUE  |
| Veillonellales                            | increased | increased | TRUE  |
| Victivallales                             | increased | decreased | FALSE |
| \[Clostridium\] leptum                    | decreased | decreased | TRUE  |
| \[Clostridium\] symbiosum                 | increased | increased | TRUE  |
| \[Ruminococcus\] lactaris                 | decreased | increased | FALSE |
| uncultured Oscillospiraceae bacterium     | increased | increased | TRUE  |

For taxa reported by both sequencing platforms, this compares the
majority-reported direction from WMS studies against the
majority-reported direction from 16S studies. High agreement supports
the idea that a taxon’s association with CRC is a real biological signal
rather than a platform-specific artifact; disagreement flags taxa where
sequencing approach may be driving the result.

## Taxon enrichment analysis

``` r

n.universe <- n_distinct(taxon.long$taxon)
n.sigs     <- length(sigs.feces)

enrichment.df <- taxon.freq %>%
  distinct(taxon, total_studies) %>%
  mutate(
    observed         = total_studies,
    expected         = n.sigs / n.universe,
    enrichment_ratio = observed / expected,
    log2_enrichment  = log2(enrichment_ratio),
    p_value = purrr::map_dbl(observed, ~ binom.test(
      x = .x, n = n.sigs, p = 1 / n.universe, alternative = "greater"
    )$p.value),
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_adj, desc(log2_enrichment))

sum(enrichment.df$p_adj < 0.05)
```

    ## [1] 257

``` r

knitr::kable(head(enrichment.df, 15))
```

| taxon | total_studies | observed | expected | enrichment_ratio | log2_enrichment | p_value | p_adj |
|:---|---:|---:|---:|---:|---:|---:|---:|
| Bacteroides | 20 | 20 | 0.3406836 | 58.70550 | 5.875424 | 0 | 0 |
| Lachnospiraceae | 20 | 20 | 0.3406836 | 58.70550 | 5.875424 | 0 | 0 |
| Fusobacterium | 18 | 18 | 0.3406836 | 52.83495 | 5.723421 | 0 | 0 |
| Ruminococcus | 18 | 18 | 0.3406836 | 52.83495 | 5.723421 | 0 | 0 |
| Clostridium | 17 | 17 | 0.3406836 | 49.89968 | 5.640959 | 0 | 0 |
| Faecalibacterium | 17 | 17 | 0.3406836 | 49.89968 | 5.640959 | 0 | 0 |
| Bacillota | 16 | 16 | 0.3406836 | 46.96440 | 5.553496 | 0 | 0 |
| Coprococcus | 16 | 16 | 0.3406836 | 46.96440 | 5.553496 | 0 | 0 |
| Prevotella | 16 | 16 | 0.3406836 | 46.96440 | 5.553496 | 0 | 0 |
| Roseburia | 16 | 16 | 0.3406836 | 46.96440 | 5.553496 | 0 | 0 |
| Lachnospira | 15 | 15 | 0.3406836 | 44.02913 | 5.460386 | 0 | 0 |
| Oscillospiraceae | 15 | 15 | 0.3406836 | 44.02913 | 5.460386 | 0 | 0 |
| Parvimonas | 15 | 15 | 0.3406836 | 44.02913 | 5.460386 | 0 | 0 |
| Peptostreptococcus | 15 | 15 | 0.3406836 | 44.02913 | 5.460386 | 0 | 0 |
| Porphyromonas | 15 | 15 | 0.3406836 | 44.02913 | 5.460386 | 0 | 0 |

A binomial test compares each taxon’s observed frequency against the
frequency expected if taxa were reported at random from the full pool of
fecal CRC taxa. Taxa with a small adjusted p-value are reported far more
often than chance would predict, independent of the core-signature
concordance threshold above.

## Co-occurrence and co-exclusion

``` r

top.taxa <- taxon.freq %>%
  distinct(taxon, total_studies) %>%
  slice_max(total_studies, n = 30) %>%
  pull(taxon)

presence.mat <- taxon.long %>%
  filter(taxon %in% top.taxa) %>%
  distinct(PMID, taxon) %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = taxon, values_from = present, values_fill = 0L) %>%
  column_to_rownames("PMID") %>%
  as.matrix()

co.occ.mat <- t(presence.mat) %*% presence.mat
diag(co.occ.mat) <- 0

co.occ.pairs <- which(upper.tri(co.occ.mat), arr.ind = TRUE) %>%
  as.data.frame() %>%
  mutate(
    taxon_a   = rownames(co.occ.mat)[row],
    taxon_b   = colnames(co.occ.mat)[col],
    n_studies = co.occ.mat[cbind(row, col)]
  ) %>%
  filter(n_studies > 0) %>%
  arrange(desc(n_studies)) %>%
  slice_max(n_studies, n = 15)

knitr::kable(co.occ.pairs %>% select(taxon_a, taxon_b, n_studies))
```

| taxon_a                 | taxon_b            | n_studies |
|:------------------------|:-------------------|----------:|
| Porphyromonas           | Fusobacterium      |        11 |
| Peptostreptococcus      | Parvimonas         |        11 |
| Roseburia               | Coprococcus        |        10 |
| Fusobacterium nucleatum | Parvimonas micra   |        10 |
| Porphyromonas           | Peptostreptococcus |         9 |
| Streptococcus           | Bacteroides        |         9 |
| Bacteroides             | Roseburia          |         9 |
| Clostridia              | Lachnospiraceae    |         9 |
| Coprococcus             | Lachnospiraceae    |         9 |
| Lachnospiraceae         | Enterobacteriaceae |         9 |
| Fusobacterium           | Fusobacteriota     |         9 |
| Bacteroides             | Coprococcus        |         8 |
| Bacteroides             | Lachnospiraceae    |         8 |
| Roseburia               | Lachnospiraceae    |         8 |
| Coprococcus             | Prevotella         |         8 |
| Coprococcus             | Ruminococcus       |         8 |
| Lachnospiraceae         | Oscillospiraceae   |         8 |
| Roseburia               | Lachnospira        |         8 |
| Ruminococcus            | Lachnospira        |         8 |
| Ruminococcus            | Clostridium        |         8 |
| Lachnospira             | Clostridium        |         8 |
| Porphyromonas           | Parvimonas         |         8 |
| Fusobacterium           | Parvimonas         |         8 |
| Roseburia               | Bacillota          |         8 |
| Lachnospiraceae         | Bacillota          |         8 |
| Roseburia               | Faecalibacterium   |         8 |
| Lachnospira             | Faecalibacterium   |         8 |

``` r

co.excl <- taxon.long %>%
  filter(taxon %in% top.taxa) %>%
  distinct(PMID, taxon, direction) %>%
  inner_join(
    taxon.long %>% filter(taxon %in% top.taxa) %>% distinct(PMID, taxon, direction),
    by = "PMID", suffix = c("_a", "_b"), relationship = "many-to-many"
  ) %>%
  filter(taxon_a < taxon_b, direction_a != direction_b) %>%
  count(taxon_a, taxon_b, name = "n_coexclusion") %>%
  filter(n_coexclusion >= 2) %>%
  arrange(desc(n_coexclusion))

nrow(co.excl)
```

    ## [1] 206

``` r

knitr::kable(head(co.excl, 15))
```

| taxon_a            | taxon_b            | n_coexclusion |
|:-------------------|:-------------------|--------------:|
| Bacteroides        | Lachnospiraceae    |             7 |
| Lachnospira        | Parvimonas         |             7 |
| Actinomyces        | Faecalibacterium   |             6 |
| Actinomyces        | Lachnospira        |             6 |
| Bacillota          | Pseudomonadota     |             6 |
| Bacteroides        | Streptococcus      |             6 |
| Enterobacteriaceae | Faecalibacterium   |             6 |
| Enterobacteriaceae | Lachnospiraceae    |             6 |
| Actinomyces        | Oscillibacter      |             5 |
| Bacteroides        | Bifidobacterium    |             5 |
| Bacteroides        | Coprococcus        |             5 |
| Bacteroides        | Enterobacteriaceae |             5 |
| Bacteroides        | Odoribacter        |             5 |
| Bacteroides        | Prevotella         |             5 |
| Clostridium        | Parvimonas         |             5 |

Co-occurrence identifies taxa that tend to be reported together within
the same study, suggestive of a shared microbial community shift.
Co-exclusion identifies taxa reported in *opposite* directions within
the same study, which may point to competitive or antagonistic
relationships between taxa in the CRC gut environment.

## Cross-disease overlap

``` r

comparison.conditions <- c(
  "Crohn's disease", "Ulcerative colitis", "Inflammatory Bowel Diseases",
  "Non-alcoholic fatty liver disease", "Obesity", "Colorectal adenoma"
)

other.diseases <- dat %>% filter(Condition %in% comparison.conditions)

disease.taxa <- purrr::map_dfr(comparison.conditions, function(cond) {
  sub.df <- other.diseases %>% filter(Condition == cond)
  if (nrow(sub.df) == 0) return(NULL)
  sigs.sub <- bugsigdbr::getSignatures(sub.df, tax.id.type = "taxname")
  tibble::tibble(taxon = unique(unlist(sigs.sub)), condition = cond)
})

crc.taxa.set <- unique(unlist(sigs.feces))

overlap.summary <- purrr::map_dfr(comparison.conditions, function(cond) {
  other.set <- disease.taxa %>% filter(condition == cond) %>% pull(taxon)
  if (length(other.set) == 0) return(NULL)
  tibble::tibble(
    condition          = cond,
    n_taxa_other       = length(other.set),
    n_shared_with_crc  = length(intersect(crc.taxa.set, other.set)),
    jaccard            = round(length(intersect(crc.taxa.set, other.set)) /
                                  length(union(crc.taxa.set, other.set)), 3)
  )
})

knitr::kable(overlap.summary)
```

| condition                         | n_taxa_other | n_shared_with_crc | jaccard |
|:----------------------------------|-------------:|------------------:|--------:|
| Crohn’s disease                   |          371 |               244 |   0.236 |
| Ulcerative colitis                |          262 |               195 |   0.200 |
| Non-alcoholic fatty liver disease |          296 |               206 |   0.207 |
| Obesity                           |          569 |               356 |   0.318 |
| Colorectal adenoma                |          233 |               195 |   0.206 |

Comparing the CRC taxon set against related diseases (IBD, obesity,
NAFLD, adenoma) gives a sense of specificity: a “core” CRC signature
that overlaps heavily with several unrelated conditions is more likely
to reflect general gut dysbiosis than a disease-specific mechanism.
**Note:** the `comparison.conditions` strings should be checked against
`unique(dat$Condition)` in your own session, since exact BugSigDB
condition labels can differ slightly from common usage.

## Geographic trends

``` r

crc.filtered %>%
  distinct(PMID, `Location of subjects`) %>%
  count(`Location of subjects`, sort = TRUE) %>%
  knitr::kable()
```

| Location of subjects                                       |   n |
|:-----------------------------------------------------------|----:|
| China                                                      |  41 |
| United States of America                                   |   5 |
| Japan                                                      |   4 |
| Italy                                                      |   3 |
| Austria                                                    |   2 |
| Spain                                                      |   2 |
| Taiwan                                                     |   2 |
| Austria,Canada,China,France,Italy,United States of America |   1 |
| Belgium                                                    |   1 |
| Canada,United States of America                            |   1 |
| France                                                     |   1 |
| India                                                      |   1 |
| Indonesia                                                  |   1 |
| Jordan                                                     |   1 |
| Mali                                                       |   1 |
| Morocco                                                    |   1 |
| Norway                                                     |   1 |
| Singapore                                                  |   1 |
| Sweden                                                     |   1 |
| United States of America,France                            |   1 |
| NA                                                         |   1 |

``` r

crc.filtered %>%
  distinct(PMID, `Location of subjects`, `Sequencing type`) %>%
  count(`Location of subjects`, `Sequencing type`) %>%
  ggplot(aes(x = reorder(`Location of subjects`, n), y = n, fill = `Sequencing type`)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "CRC Fecal Studies by Reported Location",
    x = NULL, y = "Number of studies"
  ) +
  theme_minimal(base_size = 11)
```

![](Colorectal_Cancer_Analysis_vignette_files/figure-html/geographic-trends-1.png)

Geographic distribution matters for a “core” CRC signature: if the
signature is driven almost entirely by studies from one region, it is
less safe to generalize as a universal CRC marker.

## Temporal trends

``` r

crc.filtered %>%
  distinct(PMID, Year, `Sequencing type`) %>%
  count(Year, `Sequencing type`) %>%
  group_by(`Sequencing type`) %>%
  arrange(Year) %>%
  mutate(cumulative_n = cumsum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = Year, y = cumulative_n, color = `Sequencing type`)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  labs(
    title = "Cumulative CRC Fecal Studies Over Time",
    x = "Publication year", y = "Cumulative number of studies", color = "Sequencing type"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
```

![](Colorectal_Cancer_Analysis_vignette_files/figure-html/temporal-trends-1.png)

Tracking publication year shows whether the evidence base is still
growing and whether older, low-throughput studies (mostly 16S) are
gradually being supplemented or superseded by newer WMS studies.

## Diagnostic-relevance framing

``` r

known.crc.markers <- c(
  "Fusobacterium nucleatum", "Parvimonas micra", "Peptostreptococcus anaerobius",
  "Solobacterium moorei", "Gemella morbillorum", "Porphyromonas asaccharolytica",
  "Prevotella intermedia", "Alistipes finegoldii"
)

marker.pattern <- paste(known.crc.markers, collapse = "|")

diagnostic.check <- core.sig %>%
  mutate(known_marker = str_detect(taxon, regex(marker.pattern, ignore_case = TRUE))) %>%
  select(taxon, majority_dir, n_studies, concordance_pct, known_marker)

knitr::kable(diagnostic.check)
```

| taxon | majority_dir | n_studies | concordance_pct | known_marker |
|:---|:---|---:|---:|:---|
| Lachnospiraceae | decreased | 18 | 77.8 | FALSE |
| Faecalibacterium | decreased | 16 | 87.5 | FALSE |
| Clostridium | decreased | 15 | 80.0 | FALSE |
| Lachnospira | decreased | 15 | 80.0 | FALSE |
| Clostridia | decreased | 12 | 75.0 | FALSE |
| Bifidobacterium | decreased | 11 | 81.8 | FALSE |
| Faecalibacterium prausnitzii | decreased | 9 | 100.0 | FALSE |
| Blautia | decreased | 9 | 88.9 | FALSE |
| Anaerobutyricum hallii | decreased | 8 | 87.5 | FALSE |
| Lachnospira eligens | decreased | 8 | 87.5 | FALSE |
| Subdoligranulum | decreased | 8 | 87.5 | FALSE |
| Actinomycetota | decreased | 7 | 85.7 | FALSE |
| Romboutsia | decreased | 7 | 85.7 | FALSE |
| Roseburia intestinalis | decreased | 6 | 100.0 | FALSE |
| Bifidobacteriaceae | decreased | 6 | 83.3 | FALSE |
| Bifidobacteriales | decreased | 6 | 83.3 | FALSE |
| Segatella copri | decreased | 6 | 83.3 | FALSE |
| Bifidobacterium animalis | decreased | 5 | 100.0 | FALSE |
| Bifidobacterium longum | decreased | 5 | 100.0 | FALSE |
| Staphylococcus | decreased | 5 | 80.0 | FALSE |
| Agathobacter | decreased | 4 | 100.0 | FALSE |
| Agathobacter rectalis | decreased | 4 | 100.0 | FALSE |
| Anaerostipes hadrus | decreased | 4 | 100.0 | FALSE |
| Bifidobacterium adolescentis | decreased | 4 | 100.0 | FALSE |
| Cyanobacteriota | decreased | 4 | 100.0 | FALSE |
| Fusicatenibacter | decreased | 4 | 100.0 | FALSE |
| Pseudocoprococcus catus | decreased | 4 | 100.0 | FALSE |
| Tyzzerella | decreased | 4 | 100.0 | FALSE |
| Anaerotruncus | decreased | 4 | 75.0 | FALSE |
| Bacteroides intestinalis | decreased | 4 | 75.0 | FALSE |
| Bacteroides thetaiotaomicron | decreased | 4 | 75.0 | FALSE |
| Bifidobacterium catenulatum | decreased | 4 | 75.0 | FALSE |
| Clostridiales bacterium | decreased | 4 | 75.0 | FALSE |
| Eubacterium coprostanoligenes | decreased | 4 | 75.0 | FALSE |
| Leptotrichia | decreased | 4 | 75.0 | FALSE |
| Micrococcaceae | decreased | 4 | 75.0 | FALSE |
| Neisseria | decreased | 4 | 75.0 | FALSE |
| Allocoprococcus comes | decreased | 3 | 100.0 | FALSE |
| Blautia obeum | decreased | 3 | 100.0 | FALSE |
| Coprobacillus | decreased | 3 | 100.0 | FALSE |
| Roseburia faecis | decreased | 3 | 100.0 | FALSE |
| Streptococcus salivarius | decreased | 3 | 100.0 | FALSE |
| Streptococcus thermophilus | decreased | 3 | 100.0 | FALSE |
| Fusobacterium | increased | 18 | 77.8 | FALSE |
| Parvimonas | increased | 15 | 100.0 | FALSE |
| Peptostreptococcus | increased | 15 | 100.0 | FALSE |
| Porphyromonas | increased | 15 | 80.0 | FALSE |
| Fusobacterium nucleatum | increased | 14 | 92.9 | TRUE |
| Enterobacteriaceae | increased | 14 | 78.6 | FALSE |
| Pseudomonadota | increased | 12 | 75.0 | FALSE |
| Parvimonas micra | increased | 11 | 100.0 | TRUE |
| Fusobacteriota | increased | 11 | 81.8 | FALSE |
| Parabacteroides | increased | 11 | 81.8 | FALSE |
| Gemella morbillorum | increased | 10 | 100.0 | TRUE |
| Peptostreptococcus stomatis | increased | 10 | 100.0 | FALSE |
| Solobacterium moorei | increased | 9 | 100.0 | TRUE |
| Bacteroides fragilis | increased | 9 | 88.9 | FALSE |
| Eggerthella | increased | 9 | 88.9 | FALSE |
| Lactobacillus | increased | 9 | 88.9 | FALSE |
| \[Clostridium\] symbiosum | increased | 9 | 88.9 | FALSE |
| Enterococcus | increased | 9 | 77.8 | FALSE |
| Fusobacteriia | increased | 9 | 77.8 | FALSE |
| Gemella | increased | 8 | 100.0 | FALSE |
| Akkermansia | increased | 8 | 87.5 | FALSE |
| Escherichia | increased | 8 | 87.5 | FALSE |
| Akkermansia muciniphila | increased | 8 | 75.0 | FALSE |
| Enterobacterales | increased | 8 | 75.0 | FALSE |
| Escherichia/Shigella sp. | increased | 8 | 75.0 | FALSE |
| Fusobacteriaceae | increased | 8 | 75.0 | FALSE |
| Fusobacteriales | increased | 8 | 75.0 | FALSE |
| Lactobacillales | increased | 8 | 75.0 | FALSE |
| Escherichia coli | increased | 7 | 100.0 | FALSE |
| Peptostreptococcus anaerobius | increased | 7 | 100.0 | TRUE |
| Bacilli | increased | 7 | 85.7 | FALSE |
| Porphyromonas asaccharolytica | increased | 7 | 85.7 | TRUE |
| Lactobacillaceae | increased | 6 | 100.0 | FALSE |
| Prevotella intermedia | increased | 6 | 100.0 | TRUE |
| Verrucomicrobiota | increased | 6 | 100.0 | FALSE |
| Gammaproteobacteria | increased | 6 | 83.3 | FALSE |
| Bilophila wadsworthia | increased | 5 | 100.0 | FALSE |
| Flavonifractor plautii | increased | 5 | 100.0 | FALSE |
| Hungatella hathewayi | increased | 5 | 100.0 | FALSE |
| Odoribacter splanchnicus | increased | 5 | 100.0 | FALSE |
| Streptococcus anginosus | increased | 5 | 100.0 | FALSE |
| Dialister pneumosintes | increased | 5 | 80.0 | FALSE |
| Morganella | increased | 5 | 80.0 | FALSE |
| Oscillospira | increased | 5 | 80.0 | FALSE |
| Paraprevotella | increased | 5 | 80.0 | FALSE |
| Proteus | increased | 5 | 80.0 | FALSE |
| Selenomonadales | increased | 5 | 80.0 | FALSE |
| Shigella | increased | 5 | 80.0 | FALSE |
| Alistipes onderdonkii | increased | 4 | 100.0 | FALSE |
| Bacteroidaceae | increased | 4 | 100.0 | FALSE |
| Bacteroides salyersiae | increased | 4 | 100.0 | FALSE |
| Hungatella | increased | 4 | 100.0 | FALSE |
| Porphyromonas uenonis | increased | 4 | 100.0 | FALSE |
| Verrucomicrobiales | increased | 4 | 100.0 | FALSE |
| Verrucomicrobiia | increased | 4 | 100.0 | FALSE |
| uncultured Oscillospiraceae bacterium | increased | 4 | 100.0 | FALSE |
| Acinetobacter | increased | 4 | 75.0 | FALSE |
| Anaerostipes caccae | increased | 4 | 75.0 | FALSE |
| Atopobium | increased | 4 | 75.0 | FALSE |
| Bacteroides ovatus | increased | 4 | 75.0 | FALSE |
| Christensenella | increased | 4 | 75.0 | FALSE |
| Corynebacteriaceae | increased | 4 | 75.0 | FALSE |
| Desulfovibrionaceae | increased | 4 | 75.0 | FALSE |
| Eggerthella lenta | increased | 4 | 75.0 | FALSE |
| Enterobacter | increased | 4 | 75.0 | FALSE |
| Negativicutes | increased | 4 | 75.0 | FALSE |
| Phocaeicola vulgatus | increased | 4 | 75.0 | FALSE |
| Ralstonia | increased | 4 | 75.0 | FALSE |
| Streptococcus lutetiensis | increased | 4 | 75.0 | FALSE |
| Thomasclavelia | increased | 4 | 75.0 | FALSE |
| Abiotrophia | increased | 3 | 100.0 | FALSE |
| Actinomycetaceae | increased | 3 | 100.0 | FALSE |
| Aeromonas | increased | 3 | 100.0 | FALSE |
| Alistipes finegoldii | increased | 3 | 100.0 | TRUE |
| Alloprevotella tannerae | increased | 3 | 100.0 | FALSE |
| Anaerococcus vaginalis | increased | 3 | 100.0 | FALSE |
| Bacteroides eggerthii | increased | 3 | 100.0 | FALSE |
| Bifidobacterium breve | increased | 3 | 100.0 | FALSE |
| Comamonadaceae | increased | 3 | 100.0 | FALSE |
| Desulfovibrio desulfuricans | increased | 3 | 100.0 | FALSE |
| Eisenbergiella | increased | 3 | 100.0 | FALSE |
| Klebsiella oxytoca | increased | 3 | 100.0 | FALSE |
| Lachnospiraceae bacterium 3_1_57FAA_CT1 | increased | 3 | 100.0 | FALSE |
| Methanobrevibacter smithii | increased | 3 | 100.0 | FALSE |
| Olsenella | increased | 3 | 100.0 | FALSE |
| Parabacteroides distasonis | increased | 3 | 100.0 | FALSE |
| Parabacteroides merdae | increased | 3 | 100.0 | FALSE |
| Selenomonas sputigena | increased | 3 | 100.0 | FALSE |
| Solobacterium | increased | 3 | 100.0 | FALSE |
| Streptococcus oralis | increased | 3 | 100.0 | FALSE |
| Tissierellia | increased | 3 | 100.0 | FALSE |
| Verrucomicrobiaceae | increased | 3 | 100.0 | FALSE |

``` r

cat(
  "Known literature markers recovered in the core signature:",
  sum(diagnostic.check$known_marker), "of", length(known.crc.markers), "\n"
)
```

    ## Known literature markers recovered in the core signature: 8 of 8

As a lightweight sanity check (not a formal validation), the core
consensus signature was cross-referenced against a small set of taxa
repeatedly reported as CRC-associated in the oncomicrobiome literature
(e.g. *Fusobacterium nucleatum*, *Parvimonas micra*, *Peptostreptococcus
anaerobius*). Recovering these taxa from an independent, database-driven
pipeline is reassuring; missing them does not invalidate the analysis,
but is worth noting when interpreting the core signature’s clinical
relevance. The `known.crc.markers` list is illustrative rather than
exhaustive and should be adapted to whichever markers are most relevant
to your research question.

## Monte Carlo critical N

``` r

siglengths <- lengths(crc.filtered$`MetaPhlAn taxon names`)

criticalN <- getCriticalN(
  relevant.sigs = crc.filtered$`MetaPhlAn taxon names`,
  siglengths = siglengths
)

criticalN
```

    ## $critical_n
    ## 95% 
    ##  39 
    ## 
    ## $ci_lower
    ## 2.5% 
    ##   38 
    ## 
    ## $ci_upper
    ## 97.5% 
    ##    40 
    ## 
    ## $simulated_values
    ##    [1] 32 31 31 33 32 32 34 27 36 38 34 32 29 32 27 32 28 33 34 31 36 28 37 29
    ##   [25] 41 28 31 33 29 34 31 35 35 36 28 32 35 27 29 30 33 34 39 29 44 30 34 30
    ##   [49] 44 30 36 41 33 32 32 28 35 34 35 31 35 29 32 34 31 34 38 31 37 28 33 34
    ##   [73] 30 31 35 29 35 39 37 33 41 39 32 30 38 35 34 35 29 28 33 34 34 33 37 34
    ##   [97] 42 38 30 32 35 32 32 32 35 33 33 35 30 37 33 29 28 34 27 32 36 30 32 33
    ##  [121] 27 36 34 36 29 31 35 33 35 31 34 32 36 34 40 35 32 28 36 31 31 29 30 35
    ##  [145] 31 34 33 35 33 34 35 34 29 31 37 35 39 33 32 35 36 30 46 36 37 34 36 33
    ##  [169] 36 30 35 36 28 33 33 30 29 27 39 36 30 28 30 35 39 36 30 34 34 42 34 34
    ##  [193] 28 33 36 31 29 27 32 34 33 37 36 37 35 36 34 29 34 30 32 35 31 28 31 27
    ##  [217] 34 35 31 33 33 32 34 35 33 30 32 33 35 33 31 31 32 30 28 32 36 34 34 37
    ##  [241] 29 37 35 32 31 31 34 29 39 29 29 36 32 31 29 32 29 34 41 32 36 41 39 32
    ##  [265] 34 33 32 30 32 44 28 27 36 32 34 31 35 34 29 32 33 42 33 30 35 30 27 35
    ##  [289] 37 36 31 32 32 31 34 35 30 34 29 35 34 30 37 29 29 31 34 36 36 35 35 32
    ##  [313] 35 29 33 30 29 34 31 34 31 32 27 37 35 32 36 29 42 35 28 31 31 36 29 29
    ##  [337] 31 33 30 32 34 32 35 35 30 36 34 30 31 42 34 33 33 32 30 32 30 31 33 30
    ##  [361] 30 31 33 38 30 32 29 29 27 35 29 35 40 42 29 33 34 34 34 35 34 31 39 31
    ##  [385] 34 38 29 32 34 31 28 31 29 33 35 27 41 37 43 34 31 36 30 32 29 36 31 35
    ##  [409] 30 30 35 33 29 33 35 36 30 31 34 36 33 29 29 34 33 38 31 32 33 39 32 29
    ##  [433] 33 33 31 31 31 37 32 31 30 33 42 39 39 33 33 31 35 38 33 39 28 39 32 29
    ##  [457] 32 33 31 35 36 30 33 32 42 32 36 25 39 31 30 30 30 32 29 31 33 33 32 35
    ##  [481] 31 32 32 34 46 29 37 31 31 30 28 29 34 30 34 31 29 35 33 34 42 33 36 32
    ##  [505] 34 33 34 30 30 27 31 30 33 31 34 35 37 32 33 29 35 31 29 35 36 30 34 34
    ##  [529] 31 35 35 33 35 31 30 33 28 29 33 30 32 31 36 34 29 35 34 30 29 32 30 32
    ##  [553] 33 29 30 38 27 45 29 29 31 31 29 32 37 34 33 29 31 39 27 39 37 30 33 29
    ##  [577] 30 33 33 39 41 29 33 35 28 39 32 48 29 35 36 41 33 34 35 39 35 31 33 35
    ##  [601] 36 33 40 30 32 30 31 33 33 32 33 27 35 37 34 30 31 30 32 32 35 33 29 31
    ##  [625] 30 38 36 36 30 35 34 33 36 38 27 33 37 36 29 28 30 32 33 36 33 34 29 35
    ##  [649] 33 35 33 36 36 29 28 35 34 30 32 34 30 33 30 34 31 28 34 29 35 31 33 27
    ##  [673] 33 33 35 34 30 36 28 32 34 35 28 32 32 32 39 30 35 31 32 32 34 33 30 33
    ##  [697] 31 30 32 29 33 32 34 30 34 27 32 33 31 31 36 31 29 30 38 31 33 31 32 31
    ##  [721] 34 35 34 30 28 33 36 34 30 33 33 30 30 31 31 37 35 31 35 33 34 27 35 28
    ##  [745] 28 25 34 31 26 35 29 28 34 37 35 42 38 32 35 35 30 35 28 37 32 30 27 33
    ##  [769] 31 30 31 36 31 31 44 27 35 33 31 36 29 35 30 34 40 34 33 34 35 34 33 33
    ##  [793] 30 30 34 33 35 30 31 34 36 37 34 30 33 27 40 26 28 27 27 31 28 33 27 30
    ##  [817] 30 40 33 30 34 35 26 38 32 32 40 35 31 35 34 40 35 36 27 36 30 28 36 32
    ##  [841] 33 30 32 31 28 33 37 34 38 31 33 34 38 31 36 35 35 31 28 29 33 32 37 37
    ##  [865] 33 35 36 30 30 30 30 38 31 31 34 26 37 35 36 38 29 32 28 32 28 33 28 27
    ##  [889] 38 40 32 34 35 32 41 30 42 31 35 41 38 32 35 33 35 30 33 32 37 29 38 35
    ##  [913] 35 39 30 40 32 27 38 35 37 30 36 35 37 38 33 37 34 33 33 33 38 28 36 28
    ##  [937] 37 33 29 38 33 34 30 36 31 29 32 29 30 32 32 31 35 34 34 32 37 30 33 32
    ##  [961] 38 32 32 34 32 29 30 34 33 35 34 34 29 33 30 32 30 34 36 32 34 29 36 33
    ##  [985] 30 34 36 29 39 41 31 32 39 31 35 30 36 30 32 34
    ## 
    ## $alpha
    ## [1] 0.05
    ## 
    ## $quantile_label
    ## [1] "95th percentile"

The Monte Carlo simulation estimates the critical number of independent
studies required before repeated observation of a taxon is unlikely to
occur by chance.

## Plot critical N

``` r

hist(
  criticalN$simulated_values,
  breaks = 30,
  main = "Monte Carlo Simulation of Critical N",
  xlab = "Maximum frequency under the null model"
)

abline(
  v = criticalN$critical_n,
  col = "violet",
  lwd = 2,
  lty = 2
)
```

![](Colorectal_Cancer_Analysis_vignette_files/figure-html/criticalN-plot-1.png)

The simulated distribution illustrates the study count threshold
required for statistical significance under random sampling.

## Signature similarity

``` r

allsigs <- bugsigdbr::getSignatures(
  crc.filtered,
  tax.id.type = "taxname"
)

# Keep signatures containing at least two taxa
allsigs <- allsigs[sapply(allsigs, length) > 1]

length(allsigs)
```

    ## [1] 276

``` r

# Pairwise overlaps
mydists <- BugSigDBStats::calcPairwiseOverlaps(allsigs)

dim(mydists)
```

    ## [1] 6051    8

``` r

# Jaccard similarity matrix
jmat <- BugSigDBStats::calcJaccardSimilarity(allsigs)

# Keep BugSigDB study/experiment identifiers
row_labels <- sub("^bsdb:", "PMID", sub("_.*$", "", rownames(jmat)))
col_labels <- sub("^bsdb:", "PMID", sub("_.*$", "", colnames(jmat)))
```

Microbial signatures were converted into a pairwise Jaccard similarity
matrix using functions provided by BugSigDBStats. This matrix was
subsequently used for both the heatmap and hierarchical clustering
analyses.

## Heatmap

``` r

library(grid)
library(ComplexHeatmap)

# Calculate signature lengths
siglengths <- sapply(allsigs, length)

# Short labels for display
short.labels <- paste0("S", seq_len(nrow(jmat)))

# Scale label font size down as the matrix grows, so labels stay
# legible without overlapping - but never smaller than 5pt
n.sigs.heat <- nrow(jmat)
label.fontsize <- max(5, min(8, 10 - 0.03 * n.sigs.heat))

# Top annotation
ha <- HeatmapAnnotation(
  `Signature\nLength` = anno_barplot(
    siglengths,
    gp = gpar(fill = 4:(length(siglengths) + 3)),
    annotation_name_gp = gpar(
      fontsize = 10,
      lineheight = 0.8
    )
  )
)

# Left annotation
hr <- rowAnnotation(
  `Signature\nLength` = anno_barplot(
    siglengths,
    gp = gpar(fill = 4:(length(siglengths) + 3)),
    annotation_name_gp = gpar(
      fontsize = 10,
      lineheight = 0.8
    )
  )
)

# Create heatmap
hm <- Heatmap(
  jmat,

  top_annotation = ha,
  left_annotation = hr,

  # More room reserved for labels so long PMID/study identifiers
  # aren't truncated
  row_names_max_width = unit(14, "cm"),
  column_names_max_height = unit(6, "cm"),

  row_labels = row_labels,
  column_labels = col_labels,

  row_names_gp = gpar(fontsize = label.fontsize),
  column_names_gp = gpar(fontsize = label.fontsize),
  column_names_rot = 45,

  column_title = "Heatmap of Similarity of CRC Microbial Signatures",
  column_title_gp = gpar(fontsize = 15),

  heatmap_legend_param = list(
    title = "Jaccard Similarity",
    direction = "horizontal"
  )
)

# Draw heatmap
draw(
  hm,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "right",
  padding = unit(c(2, 2, 2, 10), "mm")
)
```

![](Colorectal_Cancer_Analysis_vignette_files/figure-html/heatmap-1.png)

``` r

# Save the heatmap as a high-resolution PNG, sized to match the
# on-screen version above so nothing is re-cropped on export
png(
  filename = "CRC_heatmap.png",
  width = max(14, 0.25 * ncol(jmat) + 6),
  height = max(11, 0.25 * nrow(jmat) + 4),
  units = "in",
  res = 300
)

draw(
  hm,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "right",
  padding = unit(c(2, 2, 2, 10), "mm")
)

dev.off()
```

    ## agg_png 
    ##       2

The heatmap visualizes pairwise Jaccard similarities among
CRC-associated microbial signatures, highlighting groups of studies with
similar taxonomic compositions.

## Dendrogram

``` r

# Convert similarity to distance
hc <- hclust(as.dist(1 - jmat))

# Assign labels directly on the hclust object (in the ORIGINAL row order
# of jmat) before converting to a dendrogram. as.dendrogram()/plot() then
# handle reordering internally via hc$order.
hc$labels <- row_labels

dend <- as.dendrogram(hc)

par(mar = c(4, 1, 3, 14))  # wide right margin for the horizontal labels

plot(
  dend,
  horiz = TRUE,
  main = "Hierarchical Cluster Map of CRC Microbial Signatures",
  xlab = "Distance (1 - Jaccard similarity)",
  cex.main = 1.1
)
```

![](Colorectal_Cancer_Analysis_vignette_files/figure-html/dendrogram-1.png)

hc \<- hclust(as.dist(1 - jmat)) hc\$labels \<- row_labels

dend \<- as.dendrogram(hc)

par(mar = c(4, 1, 3, 14))

plot( dend, horiz = TRUE, main = “Hierarchical Cluster Map of CRC
Microbial Signatures”, xlab = “Distance (1 - Jaccard similarity)”,
cex.main = 1.1 ) \`\`\`

Hierarchical clustering groups microbial signatures according to
taxonomic similarity, revealing patterns of shared microbial composition
across CRC studies.

## Discussion

This exploratory analysis identified 309 fecal, human-derived colorectal
cancer microbial signatures curated in BugSigDB. Beyond the taxonomic
summary, Monte Carlo critical-N estimation, and pairwise Jaccard
similarity/heatmap/dendrogram analyses, several additional layers of
evidence were used to characterize a candidate “core” CRC-associated
signature:

- A **formal core-signature threshold** (≥3 independent studies, ≥75%
  directional concordance) narrowed the taxon list to 135 taxa with
  consistent evidence.
- **WMS vs 16S concordance** checked whether that consistency holds
  across sequencing platforms, flagging taxa whose direction may be
  platform-dependent rather than biological.
- **Binomial enrichment testing** identified taxa reported far more
  often than chance, independent of the concordance threshold.
- **Co-occurrence and co-exclusion** analyses added a network view of
  which taxa move together or in opposition within the same studies.
- **Cross-disease overlap** with IBD, obesity, NAFLD, and adenoma
  assessed how specific the CRC signature is to CRC versus general gut
  dysbiosis.
- **Geographic and temporal trends** checked whether the evidence base
  is dominated by a single region or era.
- A **diagnostic-relevance check** against known oncomicrobiome markers
  (e.g. *Fusobacterium nucleatum*, *Parvimonas micra*) provided a
  lightweight sanity check on whether the database-derived signature
  recovers established biology.

Together, these results provide a multi-angle starting point for
identifying a candidate core CRC-associated gut microbial signature and
for prioritizing taxa for further validation, rather than relying on any
single metric in isolation.
