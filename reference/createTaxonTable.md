# Create a table of most frequent taxa in a data.frame

For each of the top taxa, counts how often it was reported increased vs.
decreased and runs a binomial test on those counts (null: 50/50).

## Usage

``` r
createTaxonTable(dat, n = 10, format = "Simple")
```

## Arguments

- dat:

  data.frame produced by
  [importBugSigDB](https://rdrr.io/pkg/bugsigdbr/man/importBugSigDB.html),
  subsetted as desired

- n:

  number of taxa to return (if sig.type=="both", this is the number of
  taxa to return for each direction)

- format:

  "Simple" (default) keeps only the core summary columns; any other
  value returns the full table, with those columns moved to the front

## Value

a data.frame with one row per taxon, giving its taxonomic level,
total/increased/decreased signature counts, and binomial test p-value.
Can be wrapped in kable_styling(kbl(...)) to format nicely.

## Examples

``` r
full.dat <- bugsigdbr::importBugSigDB()
#> Using cached version from 2026-11-09 16:28:46
createTaxonTable(full.dat, n=20)
#> # A tibble: 20 × 6
#>    `Taxon Name`     `Taxonomic Level` total_signatures increased_signatures
#>    <chr>            <chr>                        <int>                <int>
#>  1 Bacteroides      genus                          597                  301
#>  2 Streptococcus    genus                          580                  376
#>  3 Lachnospiraceae  family                         579                  268
#>  4 Prevotella       genus                          568                  296
#>  5 Lactobacillus    genus                          531                  305
#>  6 Bacillota        phylum                         459                  237
#>  7 Faecalibacterium genus                          442                  161
#>  8 Oscillospiraceae family                         441                  220
#>  9 Clostridium      genus                          410                  224
#> 10 Bifidobacterium  genus                          408                  213
#> 11 Blautia          genus                          384                  180
#> 12 Pseudomonadota   phylum                         366                  238
#> 13 Bacteroidota     phylum                         365                  180
#> 14 Actinomycetota   phylum                         362                  190
#> 15 Ruminococcus     genus                          353                  155
#> 16 Roseburia        genus                          345                  118
#> 17 Eubacteriales    order                          342                  160
#> 18 Veillonella      genus                          337                  215
#> 19 Parabacteroides  genus                          319                  178
#> 20 Haemophilus      genus                          318                  152
#> # ℹ 2 more variables: decreased_signatures <int>, `Binomial Test pval` <dbl>
```
