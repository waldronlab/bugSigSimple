# Create a table of all studies currently in data.frame

Create a table of all studies currently in data.frame

## Usage

``` r
createStudyTable(bsdb.df, includeAlso = NULL)
```

## Arguments

- bsdb.df:

  `data.frame` produced by
  [importBugSigDB](https://rdrr.io/pkg/bugsigdbr/man/importBugSigDB.html),
  pre-filtered as desired

- includeAlso:

  `character` with column names to additionally include in the output
  table (default = `NULL`)

## Value

a data.frame of basic study information. Can be wrapped in
kable_styling(kbl(...)) to format nicely.

## Examples

``` r
full.dat <- bugsigdbr::importBugSigDB()
#> Using cached version from 2026-11-09 16:28:46
createStudyTable(full.dat)
#> # A tibble: 1,246 × 9
#>    `Study Identifier` MaxCases MaxControls `Study design` Condition N_signatures
#>    <chr>                 <dbl>       <dbl> <chr>          <chr>            <int>
#>  1 AartsE._2017             19          77 case-control   Attentio…            2
#>  2 AbangeW.B._2021          87          82 cross-section… HIV infe…            3
#>  3 AbduweliUyghurtur…        4           4 prospective c… Primary …            4
#>  4 AbotsiR.E._2023         174         174 randomized co… Response…            5
#>  5 AbrahamssonT.R._2…       20          20 prospective c… Atopic e…            3
#>  6 AdambergK._2018           9           9 case-control   Obesity              6
#>  7 Adu-OppongB._2022        51          51 cross-section… Urinary …            7
#>  8 AguilarT._2020           16          56 cross-section… Obesity              2
#>  9 AhmedS._2025             NA          NA cross-section… White ma…           24
#> 10 AhnJ._2013               47          94 case-control   Colorect…            1
#> # ℹ 1,236 more rows
#> # ℹ 3 more variables: PMID <chr>, DOI <chr>, URL <chr>
```
