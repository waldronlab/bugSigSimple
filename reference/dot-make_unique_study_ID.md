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
#> Using cached version from 2026-11-20 02:47:35
bugSigSimple:::.make_unique_study_ID(full.dat[1:10, ])
#> # A tibble: 10 × 58
#>    `Study Identifier` `BSDB ID`    Study     `Study design`     PMID DOI   URL  
#>    <chr>              <chr>        <chr>     <chr>             <dbl> <chr> <chr>
#>  1 Feehan_2021        bsdb:513/1/1 Study 513 cross-sectional …    NA http… http…
#>  2 Feehan_2021        bsdb:513/2/1 Study 513 cross-sectional …    NA http… http…
#>  3 Feehan_2021        bsdb:513/3/1 Study 513 cross-sectional …    NA http… http…
#>  4 NA_NA              bsdb:562/1/1 Study 562 meta-analysis        NA http… http…
#>  5 NA_NA              bsdb:562/1/3 Study 562 meta-analysis        NA http… http…
#>  6 NA_NA              bsdb:562/2/1 Study 562 meta-analysis        NA http… http…
#>  7 NA_NA              bsdb:562/2/2 Study 562 meta-analysis        NA http… http…
#>  8 NA_NA              bsdb:562/3/1 Study 562 meta-analysis        NA http… http…
#>  9 NA_NA              bsdb:562/3/2 Study 562 meta-analysis        NA http… http…
#> 10 NA_NA              bsdb:562/4/1 Study 562 meta-analysis        NA http… http…
#> # ℹ 51 more variables: `Authors list` <chr>, Title <chr>, Journal <chr>,
#> #   Year <dbl>, Keywords <chr>, `Curation method` <lgl>,
#> #   `Curation protocol` <lgl>, submittedBy <chr>, submissionDate <chr>,
#> #   `Review state` <chr>, `Quality control tag` <lgl>, Experiment <chr>,
#> #   `Location of subjects` <chr>, `Host species` <chr>, `Body site` <chr>,
#> #   `UBERON ID` <chr>, Condition <chr>, `EFO ID` <chr>, `Group 0 name` <chr>,
#> #   `Group 1 name` <chr>, `Group 1 definition` <chr>, …
```
