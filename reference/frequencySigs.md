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
#> Using cached version from 2026-11-20 02:47:35
dat.select <- bugsigdbr::getSignatures(dat, tax.level = "genus", exact.tax.level=TRUE)
dat.select <- dat.select[grep("UP", names(dat.select))] #only "UP" signatures
frequencySigs(sigs=dat.select, n = 10)
#>   1301    816    838   1578   1678  29465   1485    848 375288 572511 
#>    592    500    496    484    383    326    304    296    285    285 
```
