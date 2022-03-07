# bugSigSimple

This repo demonstrates several analyses of curated microbe signatures from
https://bugsigdb.org. See "Articles" for existing analyses, and "Reference" for 
functions provided by this packages. To install this bugSigSimple and use these
functions, do:

```r
BiocManager::install("waldronlab/bugSigSimple")
```

In particular, look at:

1. [Delivery Mode and the Meconium Microbiome](https://waldronlab.io/bugSigSimple/articles/c-section_meconium_shaimaa.html) by Shaimaa Elsafoury for a basic creation of a study table and tables of the frequency of the most commonly identified taxa.
2. [The COVID-19 associated microbiome](http://waldronlab.io/bugSigSimple/articles/capstoneanalysis_clare.html) by Clare Grieve: for professional tables, clustered heatmaps and interactive heatmaps for identifying similar signatures, and hierarchical dendrograms.
3. [An analysis of Most frequent taxa in Major Depression and Bipolar Disorder](http://waldronlab.io/bugSigSimple/articles/capstoneanalysis_fatima.html) by Fatima Azhar: This analysis does not take advantage of current functionality provided by this package and by BugSigDBStats, but it demonstrates a manual binomial test that can be used to test $H_0$: a taxon is equally probable to be reported as increased abundance or decreased abundance. This is only a sufficiently powered test when the same taxon is reported in _at least_ 6 different studies (`binom.test(0, 6) -> p=0.03125)`)
4. [The Irritable Bowel Syndrome-associated Microbiome](http://waldronlab.io/bugSigSimple/articles/capstoneanalysis_kweku.html) by Kweku Amoo: Similar to 1 in another context.
5. [The Endometriosis-associated Microbiome](http://waldronlab.io/bugSigSimple/articles/fieldworkanalysis_samara.html) by Samara Khan: This analysis generates a hypothesis based on frequently identified genera in bugsigdb.org, then tests the hypothesis with individual-participant data from curatedMetagenomicData, using linear regression, t-test, a Poisson model, and (zero-inflated) negative binomial log-linear models. 
