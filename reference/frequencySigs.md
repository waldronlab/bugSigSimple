# Identify the most frequently recurring taxa in bugsigdb

Identify the most frequently recurring taxa in bugsigdb

## Usage

``` r
frequencySigs(sigs, n = 10)
```

## Arguments

- sigs:

  a list of signatures, such as produced by
  [getSignatures](https://rdrr.io/pkg/bugsigdbr/man/getSignatures.html)

- n:

  Number of top taxa to return (default: 10)

## Value

a named table of the top `n` most recurrent taxa (taxon names as names,
counts as values), sorted in decreasing order

## Examples

``` r
dat <- bugsigdbr::importBugSigDB()
#> Using cached version from 2026-11-09 16:28:46
dat.select <- bugsigdbr::getSignatures(dat, tax.level = "genus", exact.tax.level=TRUE)
dat.select <- dat.select[grep("UP", names(dat.select))] #only "UP" signatures
frequencySigs(sigs=dat.select, n = 10)
#>   1301   1578    816    838   1485  29465   1678 572511 239934 375288 
#>    375    305    301    296    224    215    213    180    179    178 
```
