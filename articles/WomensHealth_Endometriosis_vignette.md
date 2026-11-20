# Focused Sub-Analysis of Endometriosis Microbial Signatures

## Introduction

This vignette demonstrates a focused sub-analysis of microbial
signatures associated specifically with endometriosis (and the closely
related condition adenomyosis) using BugSigDB. This is the smallest of
the four sub-analyses by sample size in the General Women’s Health
dataset (42 signatures total across “Endometriosis”, “Adenomyosis”, and
“Adenomyosis,Endometriosis”), so thresholds and statistical power should
be interpreted with appropriate caution throughout.

## Download and Prepare Data

``` r

full.dat <- bugsigdbr::importBugSigDB(cache = FALSE)
```

``` r

ind1 <- lengths(full.dat[["MetaPhlAn taxon names"]]) > 0
ind2 <- lengths(full.dat[["NCBI Taxonomy IDs"]]) > 0
dat <- full.dat[ind1 & ind2, ]
```

## Select Endometriosis-Specific Signatures

``` r

endometriosis_terms <- c(
  "Endometriosis",
  "Adenomyosis",
  "Adenomyosis,Endometriosis"
)

endo_dat <- dat %>%
  filter(Condition %in% endometriosis_terms) %>%
  distinct()

nrow(endo_dat)
```

    ## [1] 68

``` r

endo_dat %>%
  count(Condition, sort = TRUE) %>%
  knitr::kable()
```

| Condition                 |   n |
|:--------------------------|----:|
| Endometriosis             |  63 |
| Adenomyosis,Endometriosis |   3 |
| Adenomyosis               |   2 |

## Study Characteristics

``` r

study.summary <- endo_dat %>%
  count(`Body site`, `Sequencing type`)

knitr::kable(study.summary)
```

| Body site | Sequencing type | n |
|:---|:---|---:|
| Cervical cavity | 16S | 3 |
| Endometrium | 16S | 3 |
| Feces | 16S | 25 |
| Feces | WMS | 4 |
| Lower part of vagina,Upper part of vagina,Cervical mucus | 16S | 2 |
| Mouth | 16S | 5 |
| Ovary | 16S | 2 |
| Posterior fornix of vagina | 16S | 4 |
| Upper part of vagina | 16S | 2 |
| Uterine cervix | 16S | 1 |
| Uterovesical pouch | 16S | 3 |
| Uterus | 16S | 4 |
| Vagina | 16S | 8 |
| Vaginal fluid,Feces,Endocervix | 16S | 2 |

``` r

endo_dat %>%
  count(Year) %>%
  arrange(Year) %>%
  knitr::kable()
```

| Year |   n |
|-----:|----:|
| 2016 |   5 |
| 2018 |   5 |
| 2019 |   7 |
| 2020 |  15 |
| 2021 |  12 |
| 2022 |   2 |
| 2024 |  18 |
| 2025 |   4 |

``` r

endo_dat %>%
  count(`Body site`, sort = TRUE) %>%
  ggplot(aes(x = reorder(`Body site`, n), y = n)) +
  geom_col(fill = "#FFDE00", width = 0.6) +
  coord_flip() +
  labs(
    title = "Endometriosis Signatures by Body Site",
    x = NULL, y = "Number of Signatures"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
```

![](WomensHealth_Endometriosis_vignette_files/figure-html/body-site-plot-1.png)

``` r

endo_dat %>%
  distinct(PMID, `Location of subjects`, `Sequencing type`) %>%
  filter(!is.na(`Location of subjects`)) %>%
  count(`Location of subjects`, `Sequencing type`) %>%
  ggplot(aes(x = reorder(`Location of subjects`, n), y = n, fill = `Sequencing type`)) +
  geom_col(width = 0.6) +
  coord_flip() +
  # Adds your pastel dark brown palette manually
  scale_fill_manual(values = c("#8B7365", "#A38D80", "#BCABA2")) + 
  labs(
    title = "Endometriosis Studies by Reported Location",
    x = NULL, y = "Number of Independent Studies", fill = "Platform"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
```

![](WomensHealth_Endometriosis_vignette_files/figure-html/geographic-trends-1.png)

``` r

endo_dat %>%
  distinct(PMID, Year, `Sequencing type`) %>%
  filter(!is.na(Year)) %>%
  count(Year, `Sequencing type`) %>%
  group_by(`Sequencing type`) %>%
  arrange(Year) %>%
  mutate(cumulative_n = cumsum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = Year, y = cumulative_n, color = `Sequencing type`)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  labs(
    title = "Cumulative Endometriosis Studies Over Time",
    x = "Publication year", y = "Cumulative number of studies", color = "Sequencing type"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
```

![](WomensHealth_Endometriosis_vignette_files/figure-html/temporal-trends-1.png)

## Taxonomic Summary

``` r

taxa.table <- createTaxonTable(endo_dat)
knitr::kable(head(taxa.table, 20))
```

| Taxon Name | Taxonomic Level | total_signatures | increased_signatures | decreased_signatures | Binomial Test pval |
|:---|:---|---:|---:|---:|---:|
| Lactobacillus | genus | 9 | 6 | 3 | 0.510 |
| Streptococcus | genus | 7 | 6 | 1 | 0.130 |
| Prevotella | genus | 7 | 4 | 3 | 1.000 |
| Bacteroidales | order | 6 | 2 | 4 | 0.690 |
| Atopobium | genus | 5 | 2 | 3 | 1.000 |
| Coriobacteriaceae | family | 5 | 4 | 1 | 0.370 |
| Enterococcus | genus | 5 | 5 | 0 | 0.062 |
| Megasphaera | genus | 5 | 2 | 3 | 1.000 |
| Bacteroides | genus | 5 | 3 | 2 | 1.000 |
| Oscillospira | genus | 4 | 4 | 0 | 0.120 |

## Split by Sequencing Type

``` r

endo.wms <- endo_dat %>% filter(`Sequencing type` == "WMS")
endo.16s <- endo_dat %>% filter(`Sequencing type` == "16S")

cat("WMS  - signatures:", nrow(endo.wms), "| studies:", n_distinct(endo.wms$PMID), "\n")
```

    ## WMS  - signatures: 4 | studies: 1

``` r

cat("16S  - signatures:", nrow(endo.16s), "| studies:", n_distinct(endo.16s$PMID), "\n")
```

    ## 16S  - signatures: 64 | studies: 20

## Long-Format Taxon Table

``` r

sigs.endo <- bugsigdbr::getSignatures(endo_dat, tax.id.type = "taxname")
sigs.all <- bugsigdbr::getSignatures(dat, tax.id.type = "taxname")

taxon.long <- purrr::map_dfr(seq_along(sigs.all), function(i) {
  taxa <- sigs.all[[i]]
  if (length(taxa) == 0) return(NULL)
  tibble::tibble(sig_id = i, taxon = taxa)
}) %>%
  left_join(
    dat %>%
      mutate(sig_id = row_number()) %>%
      select(sig_id, PMID, Condition, `Sequencing type`, `Abundance in Group 1`),
    by = "sig_id"
  ) %>%
  rename(direction = `Abundance in Group 1`)

nrow(taxon.long)
```

    ## [1] 96257

## Taxon Frequency Across Studies

``` r

taxon.freq <- taxon.long %>%
  filter(Condition %in% endometriosis_terms) %>%
  distinct(PMID, taxon, direction) %>%
  count(taxon, direction, name = "n_studies") %>%
  group_by(taxon) %>%
  mutate(total_studies = sum(n_studies)) %>%
  ungroup() %>%
  arrange(desc(total_studies))

n_distinct(taxon.freq$taxon)
```

    ## [1] 304

``` r

taxon.freq %>%
  distinct(taxon, total_studies) %>%
  slice_max(total_studies, n = 15) %>%
  knitr::kable()
```

| taxon                        | total_studies |
|:-----------------------------|--------------:|
| Bacteroides                  |             5 |
| Clostridium                  |             4 |
| Eubacteriales                |             4 |
| Faecalibacterium             |             4 |
| Lachnospiraceae              |             4 |
| Pseudomonadota               |             4 |
| Sphingomonas                 |             4 |
| Bacillota                    |             3 |
| Faecalibacterium prausnitzii |             3 |
| Flavonifractor               |             3 |
| Fusicatenibacter             |             3 |
| Lachnoclostridium            |             3 |
| Lactobacillales              |             3 |
| Oscillospiraceae             |             3 |
| Roseburia inulinivorans      |             3 |
| Thomasclavelia               |             3 |

## Core Consensus Signature

``` r

concordance <- taxon.long %>%
  filter(Condition %in% endometriosis_terms) %>%
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

# Given the small sample size (42 signatures total), the standard >=3/>=75% threshold
# may return very few or zero rows. Loosen to >=2 studies / >=70% concordance if so
# (matching the Cervical Cancer approach), and state that explicitly below.
core.sig <- concordance %>%
  filter(n_studies >= 2, concordance_pct >= 70) %>%
  arrange(majority_dir, desc(n_studies))

nrow(core.sig)
```

    ## [1] 29

``` r

knitr::kable(core.sig)
```

| taxon | n_studies | n_increased | n_decreased | concordance_pct | majority_dir |
|:---|---:|---:|---:|---:|:---|
| Clostridium | 3 | 1 | 3 | 100 | decreased |
| Bacilli | 2 | 0 | 2 | 100 | decreased |
| Prevotella | 2 | 0 | 2 | 100 | decreased |
| Roseburia inulinivorans | 2 | 1 | 2 | 100 | decreased |
| Eubacteriales | 4 | 4 | 0 | 100 | increased |
| Lachnospiraceae | 4 | 4 | 0 | 100 | increased |
| Alistipes finegoldii | 2 | 2 | 0 | 100 | increased |
| Bergeyella | 2 | 2 | 0 | 100 | increased |
| Blautia | 2 | 2 | 0 | 100 | increased |
| Clostridioides difficile | 2 | 2 | 0 | 100 | increased |
| Desulfovibrio desulfuricans | 2 | 2 | 0 | 100 | increased |
| Dialister | 2 | 2 | 0 | 100 | increased |
| Eggerthella lenta | 2 | 2 | 0 | 100 | increased |
| Enterococcus faecalis | 2 | 2 | 0 | 100 | increased |
| Enterococcus faecium | 2 | 2 | 0 | 100 | increased |
| Ethanoligenens harbinense | 2 | 2 | 0 | 100 | increased |
| Fretibacterium fastidiosum | 2 | 2 | 0 | 100 | increased |
| Fusicatenibacter | 2 | 2 | 1 | 100 | increased |
| Lactobacillaceae | 2 | 2 | 0 | 100 | increased |
| Lactobacillus | 2 | 2 | 0 | 100 | increased |
| Mediterraneibacter torques | 2 | 2 | 0 | 100 | increased |
| Megasphaera | 2 | 2 | 0 | 100 | increased |
| Neisseriaceae | 2 | 2 | 0 | 100 | increased |
| Oscillibacter valericigenes | 2 | 2 | 0 | 100 | increased |
| Oscillospiraceae | 2 | 2 | 1 | 100 | increased |
| Romboutsia | 2 | 2 | 0 | 100 | increased |
| Slackia heliotrinireducens | 2 | 2 | 0 | 100 | increased |
| Subdoligranulum | 2 | 2 | 0 | 100 | increased |
| Veillonella | 2 | 2 | 0 | 100 | increased |

## Core Signature Visualisation

``` r

if (nrow(core.sig) > 0) {
  core.sig %>%
    group_by(majority_dir) %>%
    slice_max(n_studies, n = 15, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      direction_label = if_else(majority_dir == "increased", "Increased", "Decreased"),
      taxon_label = paste0(str_wrap(taxon, width = 28), " (", n_studies, ")")
    ) %>%
    ggplot(aes(x = reorder(taxon_label, n_studies), y = n_studies, fill = direction_label)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = paste0(concordance_pct, "%")), hjust = -0.15, size = 3.3) +
    coord_flip(clip = "off") +
    facet_wrap(~direction_label, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Increased" = "#FFF2A3", "Decreased" = "#4A2E1B")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(
      title = "Core Endometriosis Consensus Signature",
      subtitle = "Taxa in >=2 studies with >=70% direction concordance",
      x = NULL, y = "Number of independent studies", fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
} else {
  cat("No taxa met the core-signature threshold at this sample size.")
}
```

![](WomensHealth_Endometriosis_vignette_files/figure-html/core-signature-plot-1.png)

## Taxon Enrichment Analysis

``` r

n.universe <- n_distinct(taxon.long %>% filter(Condition %in% endometriosis_terms) %>% pull(taxon))
n.sigs     <- length(sigs.endo)

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

    ## [1] 16

``` r

knitr::kable(head(enrichment.df, 15))
```

| taxon | total_studies | observed | expected | enrichment_ratio | log2_enrichment | p_value | p_adj |
|:---|---:|---:|---:|---:|---:|---:|---:|
| Bacteroides | 5 | 5 | 0.2105263 | 23.75 | 4.569856 | 0.0000025 | 0.0007595 |
| Clostridium | 4 | 4 | 0.2105263 | 19.00 | 4.247928 | 0.0000635 | 0.0027598 |
| Eubacteriales | 4 | 4 | 0.2105263 | 19.00 | 4.247928 | 0.0000635 | 0.0027598 |
| Faecalibacterium | 4 | 4 | 0.2105263 | 19.00 | 4.247928 | 0.0000635 | 0.0027598 |
| Lachnospiraceae | 4 | 4 | 0.2105263 | 19.00 | 4.247928 | 0.0000635 | 0.0027598 |
| Pseudomonadota | 4 | 4 | 0.2105263 | 19.00 | 4.247928 | 0.0000635 | 0.0027598 |
| Sphingomonas | 4 | 4 | 0.2105263 | 19.00 | 4.247928 | 0.0000635 | 0.0027598 |
| Bacillota | 3 | 3 | 0.2105263 | 14.25 | 3.832890 | 0.0012765 | 0.0242539 |
| Faecalibacterium prausnitzii | 3 | 3 | 0.2105263 | 14.25 | 3.832890 | 0.0012765 | 0.0242539 |
| Flavonifractor | 3 | 3 | 0.2105263 | 14.25 | 3.832890 | 0.0012765 | 0.0242539 |
| Fusicatenibacter | 3 | 3 | 0.2105263 | 14.25 | 3.832890 | 0.0012765 | 0.0242539 |
| Lachnoclostridium | 3 | 3 | 0.2105263 | 14.25 | 3.832890 | 0.0012765 | 0.0242539 |
| Lactobacillales | 3 | 3 | 0.2105263 | 14.25 | 3.832890 | 0.0012765 | 0.0242539 |
| Oscillospiraceae | 3 | 3 | 0.2105263 | 14.25 | 3.832890 | 0.0012765 | 0.0242539 |
| Roseburia inulinivorans | 3 | 3 | 0.2105263 | 14.25 | 3.832890 | 0.0012765 | 0.0242539 |

## Monte Carlo Critical N

``` r

siglengths <- lengths(endo_dat$`MetaPhlAn taxon names`)

criticalN <- getCriticalN(
  relevant.sigs = endo_dat$`MetaPhlAn taxon names`,
  siglengths = siglengths
)

criticalN
```

    ## $critical_n
    ## 95% 
    ##  13 
    ## 
    ## $ci_lower
    ## 2.5% 
    ##   13 
    ## 
    ## $ci_upper
    ## 97.5% 
    ##    13 
    ## 
    ## $simulated_values
    ##    [1] 12  9 10 10  7 11  8  7  9  8  9  8  7 10  7 10 14  8  9 12  9  9  8  9
    ##   [25]  9 12  9 12  8  9  7  8 11  9  9  9  9 11  9 10  8 10 11  8 10 13 10 11
    ##   [49] 10 10 11 11 10  8 10 11 11 13 12  8 10  9  9  8 10  9 12 13 12 11 10  9
    ##   [73] 13  8  8 12  9 10  7 11 10 12  9 11  9  8  9 10 11 10  9  9 11  8 11 11
    ##   [97]  9 13 13  9  9  9  8  9  7 10 13  9  8  9  7  9  9 10  6  8  8 14  8  9
    ##  [121]  9 11 10  9 11  9  9 10  8 10 11  8 12  9  8 11  8  8  9  9 10 14 10 12
    ##  [145] 12 10 10 10  9 10  6 10 11  9  9  9  8 11  9 12 11  8  9 14  9 10 10 11
    ##  [169] 11  7  9  9 16 11  9  9 10 13 10  8  8 12 14  7 11 14 12  8  7 10 11 15
    ##  [193]  8  8  9 11 13  8 10 10 12 11 12  9  6 10  8 12  8 13 10 11  9  8  9 10
    ##  [217]  9  9 11  9 10 12  8 12 10  9 10 10 11  9 10  8 10 12  9  9  7  9 11 10
    ##  [241] 14 11  9 10 11  9 10  8 10 10 10 13 11 11  8  9  8 10  9 13  9 10 13  7
    ##  [265]  8  9 10  9  7  8  8  9  8  9  9  8 12 10 15  9 12  9  9 10  9 10 12  9
    ##  [289] 14  8  8  9  7  8 10  8 16 10 10  9 12  8 10  9  8 10  9  8  9  9 10  8
    ##  [313]  9 10 10 11 13 11 11  8  9  9 10 10  9  9 11 12 10 10 10 13 10 11 10 12
    ##  [337] 11  7 11  9 12 10 11 12 13 13 11  9 13 10  9  7  9 10  9 10 10  7 11  9
    ##  [361] 11  7 12 10  9  9  9 10  9 11 11 10 11 10 12 12  8  9 14 10 10 13 14  9
    ##  [385] 10 12  8  7 11  9 10  7 13 10 11  8  9  9 10  9 12 13  8 10  8 12  8 11
    ##  [409]  8 11 14  8  9  9 13 10 15 11 10 11  8  8 10 11 14  9  8  9  9 12  9 11
    ##  [433]  9 11  9  9 11 10  9  8 13 11  8 16 11 11 11  7  8 10 10  8 16 12  9 10
    ##  [457]  9 10  8 11 10  8  8  9 13 10 14  9  9  9 11  8 12 10  9 12  9 10  9 10
    ##  [481] 13  8  9 13  9 11 12 10 12 15 10  8 10 11  7  8 10  9  8  9 10  8 10  7
    ##  [505]  7 16 12  8  9 10  8 11  9 10 12 11  8 10  9  8 10  8  9 10  8  7 10 10
    ##  [529] 11  8 11  9  9 11 10  9  8 10  9 11 14  8  9  9 15  8  8 13  8 13  7 11
    ##  [553] 12  8  9  9 12 13 10  9 10 10 12 13 10 10  9 10  9  8  9  9  9 10  8 10
    ##  [577] 14 10 10  8 12 10 10  8 13  9 10  9 11  9  9 12 10 13 11  8  9  9 10 10
    ##  [601]  9 10  8  8  8  9  9  8  9 12 10 11 10  8 16  9 12 11  9  8  9 11  8 13
    ##  [625]  8 10  9  9  9 10 10 12  9  8  8  9  8 11 11  9  8  9  9 10  8 11  7 10
    ##  [649]  8  9 10  9  9  7 10 10 12  9  9  9 12 11 10 11 12 13  9  8  9 14  8 10
    ##  [673] 12  9 10  8 12 10  9 10 10  9 12  8 11 10 17  7 10 10 11  9 10  9 12  8
    ##  [697]  8 12 10  9 12  8 10  9 11  8  8 11 12 11  6  8  9  9 11 11 10 11 10  9
    ##  [721] 10  9  8 10  9 12  9  9 12  9  9  9 14 10 12  9  8 11 12 10 11 10  9 12
    ##  [745] 11  9 12 11  9 10  7  8  9 10  9  8 11  9 11 12 11 12 10  9  9  8 10 12
    ##  [769]  8  9 11 11 10 10  8 11  9  8 10  8  6 12 11  8  9 10 11  9 10  9  7 10
    ##  [793] 12  9  9 10 13 10  8  7 12  9  9  9  8 10  8  9 10 12  8 10  8  7  8  9
    ##  [817]  9  9 10 10 10 10 12 11 10 11  8 13  9  9  9  9 10 12 11  9 12 12  8 13
    ##  [841]  8 11 12 11  9 11 11  9 10  9 13 10 11  9  7  9 13  9 12  9 11 12 11 10
    ##  [865] 11  8  8 10  9 11  8 13  9  9 14  9 10  9 19  9 11  8  9 11  9 12 11 10
    ##  [889]  8  9  9  9  9 11 10  9  9 13 10  9  9 10  9 10 13 13 10 12 10  9 11  9
    ##  [913]  8 11  9 11  7  9 10  8 11 14  8  9 12 11  7  8  9  9 10 10  9  9 11 10
    ##  [937]  8  8  9 11  7 11  9  9  9 12 10 11 13 10  9 12 10 12 14  9 11  8 10 11
    ##  [961] 17 10  9  9  9 11 12 10 10 12 11 10 10  9 12  8  9 10  8 12  9 10  9  8
    ##  [985]  7  9  9 10 10 10  9  9  9  9  9  8 10 10  7  9
    ## 
    ## $alpha
    ## [1] 0.05
    ## 
    ## $quantile_label
    ## [1] "95th percentile"

``` r

hist(
  criticalN$simulated_values,
  breaks = 30,
  main = "Monte Carlo Simulation of Critical N (Endometriosis)",
  xlab = "Maximum frequency under the null model"
)

abline(v = criticalN$critical_n, col = "violet", lwd = 2, lty = 2)
```

![](WomensHealth_Endometriosis_vignette_files/figure-html/criticalN-plot-1.png)

## Cross-Condition Taxa Overlap Analysis

``` r

endo_taxa <- taxon.long %>%
  filter(Condition %in% endometriosis_terms) %>%
  distinct(taxon) %>%
  pull(taxon)

shared_taxa_counts <- taxon.long %>%
  filter(taxon %in% endo_taxa) %>%
  group_by(taxon) %>%
  summarise(
    total_womens_health_studies = n_distinct(PMID),
    distinct_conditions_linked = n_distinct(Condition),
    top_other_conditions = paste(
      head(setdiff(unique(Condition), endometriosis_terms), 3),
      collapse = ", "
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(distinct_conditions_linked))

shared_taxa_counts %>%
  slice_max(distinct_conditions_linked, n = 20, with_ties = FALSE) %>%
  knitr::kable(
    caption = "Top Endometriosis-Associated Taxa and Their Links to Other Conditions",
    col.names = c("Taxon Name", "Total WH Studies", "Total Distinct Conditions", "Other Conditions it Often Appears In")
  )
```

| Taxon Name | Total WH Studies | Total Distinct Conditions | Other Conditions it Often Appears In |
|:---|---:|---:|:---|
| Streptococcus | 557 | 334 | COVID-19, Health study participation, Diet |
| Bacteroides | 548 | 321 | Health study participation, Irritable bowel syndrome, Autism |
| Prevotella | 542 | 320 | COVID-19, Health study participation, Illness severity status |
| Lachnospiraceae | 508 | 302 | Type II diabetes mellitus, Autoimmune disease, Smoking behavior |
| Oscillospiraceae | 447 | 300 | Autism, Autoimmune disease, Obesity |
| Bacillota | 468 | 296 | Parkinson’s disease, Irritable bowel syndrome, Type II diabetes mellitus |
| Faecalibacterium | 411 | 287 | Health study participation, Autism, Autoimmune disease |
| Lactobacillus | 485 | 284 | Health study participation, Antimicrobial agent, Oxygen |
| Bifidobacterium | 428 | 270 | Health study participation, Parkinson’s disease, Food allergy |
| Bacteroidota | 393 | 246 | Irritable bowel syndrome, Autism, Autoimmune disease |
| Pseudomonadota | 379 | 244 | Irritable bowel syndrome, Type II diabetes mellitus, Autism |
| Blautia | 395 | 243 | Health study participation, Autoimmune disease, Irritable bowel syndrome |
| Clostridium | 371 | 238 | Health study participation, Irritable bowel syndrome, Type II diabetes mellitus |
| Veillonella | 352 | 229 | COVID-19, Health study participation, Antimicrobial agent |
| Roseburia | 347 | 221 | Health study participation, Autoimmune disease, Smoking behavior |
| Ruminococcus | 346 | 219 | Health study participation, Type II diabetes mellitus, Autoimmune disease |
| Haemophilus | 311 | 209 | Health study participation, Illness severity status, Diet |
| Fusobacterium | 306 | 207 | Health study participation, Antimicrobial agent, Autism |
| Enterobacteriaceae | 312 | 206 | Type II diabetes mellitus, Autism, Smoking behavior |
| Eubacteriales | 318 | 203 | Type II diabetes mellitus, Autism, Autoimmune disease |

Top Endometriosis-Associated Taxa and Their Links to Other Conditions
{.table}

## Signature Similarity and Heatmap

``` r

allsigs <- bugsigdbr::getSignatures(endo_dat, tax.id.type = "taxname")
allsigs <- allsigs[sapply(allsigs, length) > 1]

if (length(allsigs) > 1) {
  jmat <- BugSigDBStats::calcJaccardSimilarity(allsigs)
  siglengths <- sapply(allsigs, length)

  ha <- HeatmapAnnotation(
    `Signature\nLength` = anno_barplot(siglengths, gp = gpar(fill = "#4682B4"), height = unit(1.2, "cm"))
  )
  hr <- rowAnnotation(
    `Signature\nLength` = anno_barplot(siglengths, gp = gpar(fill = "#4682B4"), width = unit(1.2, "cm"))
  )

  hm <- Heatmap(
    jmat,
    top_annotation = ha,
    left_annotation = hr,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 45,
    column_title = "Heatmap of Endometriosis Microbial Signatures",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(title = "Jaccard Similarity", direction = "horizontal", legend_width = unit(6, "cm"))
  )

  draw(hm, heatmap_legend_side = "bottom", padding = unit(c(5, 5, 5, 10), "mm"))

  png(filename = "Endometriosis_heatmap.png", width = 10, height = 9, units = "in", res = 300)
  draw(hm, heatmap_legend_side = "bottom", padding = unit(c(5, 5, 5, 10), "mm"))
  dev.off()
} else {
  cat("Fewer than 2 multi-taxon signatures available - heatmap not generated.")
}
```

![](WomensHealth_Endometriosis_vignette_files/figure-html/heatmap-1.png)

    ## agg_png 
    ##       2

Note: with only 42 total signatures for this sub-analysis, the heatmap
is small enough that row/column labels are kept visible (unlike
General/Breast Cancer, which had to hide labels given hundreds of
signatures) - adjust `row_names_gp`/`column_names_gp` font size further
if labels still overlap once you see the actual output.

## Dendrogram

``` r

if (exists("jmat") && nrow(jmat) > 1) {
  hc <- hclust(as.dist(1 - jmat))
  dend <- as.dendrogram(hc)

  par(mar = c(5, 2, 4, 10))

  plot(
    dend,
    horiz = TRUE,
    main = "Hierarchical Cluster Map of Endometriosis Microbial Signatures",
    xlab = "Distance (1 - Jaccard similarity)",
    cex.main = 1.2
  )
}
```

![](WomensHealth_Endometriosis_vignette_files/figure-html/dendrogram-1.png)
