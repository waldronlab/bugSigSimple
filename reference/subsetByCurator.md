# Subset a data.frame of signatures by curator

Subset a data.frame of signatures by curator

## Usage

``` r
subsetByCurator(dat, curator, curator.column = "Curator")
```

## Arguments

- dat:

  data.frame produced by
  [importBugSigDB](https://rdrr.io/pkg/bugsigdbr/man/importBugSigDB.html)

- curator:

  character vector of one or more curator names to subset by

- curator.column:

  name of column of curators in data.frame (default: "Curator")

## Value

data.frame subsetted by curator

## Examples

``` r
full.dat <- bugsigdbr::importBugSigDB()
#> Using cached version from 2026-11-09 16:28:46
fatima.dat <- subsetByCurator(full.dat, curator="Fatima Zohra")
```
