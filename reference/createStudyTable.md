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
#> Using cached version from 2026-11-20 02:47:35
createStudyTable(full.dat)
#> # A tibble: 1,846 × 9
#>    `Study Identifier` MaxCases MaxControls `Study design` Condition N_signatures
#>    <chr>                 <dbl>       <dbl> <chr>          <chr>            <int>
#>  1 AartiDarra_2023          54          54 cross-section… Hypergly…            8
#>  2 AartsE._2017             19          77 case-control   Attentio…            2
#>  3 AbangeW.B._2021          87          82 cross-section… HIV infe…            3
#>  4 AbbondioM._2025           6           6 laboratory ex… Response…            9
#>  5 AbduweliUyghurtur…        4           4 prospective c… Primary …            4
#>  6 AbiolaIsawumi_2023       57          50 case-control   Gastroen…            4
#>  7 AbotsiR.E._2023         174         174 randomized co… Response…            5
#>  8 AbrahamssonT.R._2…       20          20 prospective c… Atopic e…            3
#>  9 AdambergK._2018           9           9 case-control   Obesity              6
#> 10 AdebayoA.S._2017         46          25 case-control   Urinary …           14
#> # ℹ 1,836 more rows
#> # ℹ 3 more variables: PMID <chr>, DOI <chr>, URL <chr>
```
