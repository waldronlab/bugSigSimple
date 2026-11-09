# Build unique study identifiers for a table of signatures

Generates a study identifier the way curatedMetagenomicData does it: the
first author from the `Authors list` column (usually surname + initials,
e.g. "ZellerG") with whitespace stripped, plus the year, e.g.
"ZellerG_2014". Duplicates get a ".1", ".2", ... suffix. Also normalizes
bare DOIs (starting with "10.") to full <https://doi.org/> links.

## Usage

``` r
.make_unique_study_ID(bsdb.df)
```

## Arguments

- bsdb.df:

  `data.frame` produced by
  [importBugSigDB](https://rdrr.io/pkg/bugsigdbr/man/importBugSigDB.html),
  pre-filtered as desired

## Value

`bsdb.df` with a new leading `Study Identifier` column

## Note

Author: Giacomo Antonello (2025-03-17)

## Examples

``` r
full.dat <- bugsigdbr::importBugSigDB()
#> Using cached version from 2026-11-09 16:28:46
bugSigSimple:::.make_unique_study_ID(full.dat[1:10, ])
#> # A tibble: 10 × 52
#>    `Study Identifier`     `BSDB ID`    Study    `Study design`  PMID DOI   URL  
#>    <chr>                  <chr>        <chr>    <chr>          <dbl> <chr> <chr>
#>  1 NA_NA                  bsdb:83/1/1  Study 83 time series /…    NA NA    NA   
#>  2 NA_NA                  bsdb:83/1/2  Study 83 time series /…    NA NA    NA   
#>  3 ChristinaPavloudi_2023 bsdb:817/1/1 Study 8… case-control      NA NA    http…
#>  4 ChristinaPavloudi_2023 bsdb:817/1/2 Study 8… case-control      NA NA    http…
#>  5 ChristinaPavloudi_2023 bsdb:817/2/1 Study 8… case-control      NA NA    http…
#>  6 ChristinaPavloudi_2023 bsdb:817/2/2 Study 8… case-control      NA NA    http…
#>  7 ChristinaPavloudi_2023 bsdb:817/3/1 Study 8… case-control      NA NA    http…
#>  8 ChristinaPavloudi_2023 bsdb:817/3/2 Study 8… case-control      NA NA    http…
#>  9 AlessandraLunis_2023   bsdb:879/1/1 Study 8… prospective c…    NA http… http…
#> 10 AlessandraLunis_2023   bsdb:879/2/1 Study 8… prospective c…    NA http… http…
#> # ℹ 45 more variables: `Authors list` <chr>, Title <chr>, Journal <chr>,
#> #   Year <dbl>, Keywords <chr>, Experiment <chr>, `Location of subjects` <chr>,
#> #   `Host species` <chr>, `Body site` <chr>, `UBERON ID` <chr>,
#> #   Condition <chr>, `EFO ID` <chr>, `Group 0 name` <chr>,
#> #   `Group 1 name` <chr>, `Group 1 definition` <chr>,
#> #   `Group 0 sample size` <dbl>, `Group 1 sample size` <dbl>,
#> #   `Antibiotics exclusion` <chr>, `Sequencing type` <chr>, …
```
