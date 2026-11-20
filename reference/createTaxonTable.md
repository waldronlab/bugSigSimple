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
#> Using cached version from 2026-11-20 02:47:35
createTaxonTable(full.dat, n=20)
#> # A tibble: 20 × 6
#>    `Taxon Name`     `Taxonomic Level` total_signatures increased_signatures
#>    <chr>            <chr>                        <int>                <int>
#>  1 Bacteroides      genus                          987                  500
#>  2 Streptococcus    genus                          945                  593
#>  3 Prevotella       genus                          926                  496
#>  4 Lachnospiraceae  family                         856                  402
#>  5 Lactobacillus    genus                          810                  484
#>  6 Bacillota        phylum                         751                  362
#>  7 Oscillospiraceae family                         736                  365
#>  8 Bifidobacterium  genus                          694                  383
#>  9 Faecalibacterium genus                          665                  241
#> 10 Blautia          genus                          623                  285
#> 11 Bacteroidota     phylum                         612                  300
#> 12 Pseudomonadota   phylum                         598                  389
#> 13 Clostridium      genus                          586                  304
#> 14 Actinomycetota   phylum                         577                  304
#> 15 Veillonella      genus                          548                  326
#> 16 Ruminococcus     genus                          539                  237
#> 17 Roseburia        genus                          533                  186
#> 18 Eubacteriales    order                          518                  245
#> 19 Parabacteroides  genus                          502                  285
#> 20 Fusobacterium    genus                          485                  296
#> # ℹ 2 more variables: decreased_signatures <int>, `Binomial Test pval` <dbl>
```
