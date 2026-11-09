# getCriticalN performs a Monte Carlo simulation to estimate the number of times the most frequent taxon is expected to be observed in a list of signatures

getCriticalN performs a Monte Carlo simulation to estimate the number of
times the most frequent taxon is expected to be observed in a list of
signatures

## Usage

``` r
getCriticalN(
  relevant.sigs,
  siglengths,
  alpha = 0.05,
  nsim = 1000,
  ci = 0.95,
  ci_nsim = 500
)
```

## Arguments

- relevant.sigs:

  a list of signatures that form the "background" from which taxa for
  simulated signatures will be drawn. These are used to estimate how
  frequently taxa occur

- siglengths:

  The sizes of signatures found in a set of related studies. Simulated
  signatures will match these in number and size.

- alpha:

  Probability at which a critical threshold will be calculated (default:
  0.05)

- nsim:

  Number of simulations (default: 1000)

- ci:

  Width of the bootstrap confidence interval to compute around the
  critical threshold, e.g. 0.95 for a 95% CI (default: 0.95). Set to
  NULL to skip.

- ci_nsim:

  Number of bootstrap resamples used to compute the CI (default: 500)

## Value

A list: `critical_n` is the 1 - alpha quantile of simulated values for
the maximum number of times any taxon is identified,
`ci_lower`/`ci_upper` give the bootstrap CI around it (NA if
`ci = NULL`), `simulated_values` holds the raw simulated maxima, plus
`alpha` and a `quantile_label` for printing.

## Details

E.g. for alpha = 0.05, we expect only a 5% chance that any taxon will be
identified N times or more.

## Examples

``` r
full.dat <- bugsigdbr::importBugSigDB()
#> Using cached version from 2026-11-09 16:28:46
my.dat <- full.dat[full.dat$Curator == "Mst Afroza Parvin", ]
relevant.dat <- full.dat[full.dat$`Body site` %in% my.dat$`Body site`, ]
relevant.sigs <- bugsigdbr::getSignatures(my.dat)
my.sigs.increased <- relevant.sigs[grep("UP", names(relevant.sigs))]
(my.siglengths <- sapply(my.sigs.increased, length))
#>                                                                                                                  bsdb:22546742/1/1_Obesity:overweight/obese_vs_control-(BMI-at-normal-range)_UP 
#>                                                                                                                                                                                               1 
#>                                                                                                                               bsdb:23459324/4/1_Obesity:overweight-subjects_vs_lean-subjects_UP 
#>                                                                                                                                                                                               1 
#>                                                                                            bsdb:23631345/1/1_Obesity:obese-(obese/overweight/morbidly-obese)_vs_control-(lean/normal-weight)_UP 
#>                                                                                                                                                                                               1 
#>                                                              bsdb:25444008/1/1_Premature-birth:2-day-old-Pre-term-infants_vs_2-day-old-Full-term-infants(born-after-uncomplicated-pregnancy)_UP 
#>                                                                                                                                                                                               7 
#>                                                                    bsdb:25444008/2/1_Premature-birth:10-day-old-Pre-term_vs_10-day-old-Full-term-infants(born-after-uncomplicated-pregnancy)_UP 
#>                                                                                                                                                                                               3 
#>                                                                                  bsdb:25444008/3/1_Premature-birth:30-and-90-day-old-Pre-term-infants_vs_30-and-90-day-old-Full-term-infants_UP 
#>                                                                                                                                                                                               1 
#>      bsdb:25444008/9/1_Antimicrobial-agent:30-day-old-pre-tem-neonates-with-mothers-exposed-to-antibiotics-(IAP)_vs_30-day-old-full-term-neonates-not-exposed-to-antibiotics-(IAP/postnatal)_UP 
#>                                                                                                                                                                                               1 
#>                                                                                                                                    bsdb:26551842/1/1_Obesity:obese-children_vs_lean-children_UP 
#>                                                                                                                                                                                               2 
#>                                                                                                                               bsdb:26551842/2/1_Obesity:overweight-children_vs_lean-children_UP 
#>                                                                                                                                                                                               2 
#>                                                                                                                                    bsdb:26551842/3/1_Obesity:obese-children_vs_lean-children_UP 
#>                                                                                                                                                                                               1 
#>                                                                                                                               bsdb:26551842/4/1_Obesity:overweight-children_vs_lean-children_UP 
#>                                                                                                                                                                                               1 
#>                                                                                                                                             bsdb:26745497/1/1_Obesity:obese_vs_normal-weight_UP 
#>                                                                                                                                                                                               2 
#>                                         bsdb:27411009/1/1_Antimicrobial-agent:oral-intake-of-vancomycin-(directed-against-Gram-positive-bacteria)_vs_or-placebo-(microcrystalline-cellulose)_UP 
#>                                                                                                                                                                                              18 
#>                                                                               bsdb:28789705/1/1_Antimicrobial-agent:Infants-to-mothers-received-IAP_vs_Infants-to-mothers-didn't-receive-IAP_UP 
#>                                                                                                                                                                                               5 
#>                                                         bsdb:28789705/2/1_Antimicrobial-agent:Breast-fed-Infants-to-mothers-received-IAP_vs_breast-fed-Infants-to-mothers-didn't-receive-IAP_UP 
#>                                                                                                                                                                                               2 
#>                                                       bsdb:28789705/3/1_Antimicrobial-agent:Formula-fed-Infants-to-mothers-received-IAP_vs_Formula-fed-Infants-to-mothers-didn't-receive-IAP_UP 
#>                                                                                                                                                                                               1 
#>                                            bsdb:28808302/1/1_Antimicrobial-agent:combination-of-penicillin-and-moxalactam-group-(postnatal-day-3)_vs_Antibiotic-free-group-(postnatal-day-3)_UP 
#>                                                                                                                                                                                               5 
#>                                                             bsdb:28808302/2/1_Antimicrobial-agent:piperacillin-tazobactam-group-(postnatal-day-3)_vs_Antibiotic-free-group-(postnatal-day-3)_UP 
#>                                                                                                                                                                                               2 
#>                                            bsdb:28808302/3/1_Antimicrobial-agent:combination-of-penicillin-and-moxalactam-group-(postnatal-day-7)_vs_Antibiotic-free-group-(postnatal-day-7)_UP 
#>                                                                                                                                                                                               2 
#>                                                             bsdb:28808302/4/1_Antimicrobial-agent:piperacillin-tazobactam-group-(postnatal-day-7)_vs_Antibiotic-free-group-(postnatal-day-7)_UP 
#>                                                                                                                                                                                               2 
#>                                                                                                                      bsdb:29935270/1/1_Obesity:overweight-children_vs_normal-weight-children_UP 
#>                                                                                                                                                                                               6 
#>                                                                                                                                      bsdb:29935270/2/1_Obesity:Overweight-children_vs_Adults_UP 
#>                                                                                                                                                                                              18 
#>                                                                                                                                         bsdb:30267022/1/1_Obesity:obese-children_vs_controls_UP 
#>                                                                                                                                                                                              10 
#>                     bsdb:31006170/1/1_Antimicrobial-agent:infants(6-weeks)-to-mothers-exposed-to-intrapartum-antibiotic_vs_infants-(6-weeks)-to-mothers-unexposed-to--intrapartum-antibiotic_UP 
#>                                                                                                                                                                                               4 
#>                       bsdb:31006170/2/1_Antimicrobial-agent:infants(1-year)-to-mothers-exposed-to-intrapartum-antibiotic_vs_infants-(1-year)-to-mothers-unexposed-to--intrapartum-antibiotic_UP 
#>                                                                                                                                                                                               4 
#>                                     bsdb:31006170/4/1_Antimicrobial-agent:infants(1-year)-to-mothers-exposed-to-the-penicillin-group_vs_infants-(1-year)-to-mothers-unexposed-to--antibiotic_UP 
#>                                                                                                                                                                                               2 
#>                                bsdb:31006170/5/1_Antimicrobial-agent:infants(6-weeks)-to-mothers-exposed-to-multi-class-antibiotics_vs_infants-(6-weeks)-to-mothers-unexposed-to--antibiotic_UP 
#>                                                                                                                                                                                               1 
#>       bsdb:31226994/1/1_Antimicrobial-agent:3-month-old-infants-exposed-to-antibiotics-in-the-second-trimester_vs_3-month-old-infants-not-exposed-to-antibiotics-during-the-second-trimester_UP 
#>                                                                                                                                                                                               2 
#> bsdb:31226994/3/1_Antimicrobial-agent:12-month-old-infants-exposed-to-antibiotics-during-the-second-trimester_vs_12-month-old-infants-not-exposed-to-antibiotics-during-the-second-trimester_UP 
#>                                                                                                                                                                                               4 
#>                                                                                         bsdb:31337807/1/1_Antimicrobial-agent:Postnatal-antibiotic-group-(2-day-old-infant)_vs_control-group_UP 
#>                                                                                                                                                                                               2 
#>                                                                                                        bsdb:31337807/3/1_Antimicrobial-agent:IAP-group-(6-month-old-infant)_vs_control-group_UP 
#>                                                                                                                                                                                               2 
#>                                                                                     bsdb:31337807/4/1_Antimicrobial-agent:IAP+postnatal-antibiotic-group-(2-day-old-infant)_vs_control-group_UP 
#>                                                                                                                                                                                               1 
#>                                                                                   bsdb:31337807/6/1_Antimicrobial-agent:IAP+postnatal-antibiotic-group-(6-month-old-infant)_vs_control-group_UP 
#>                                                                                                                                                                                               1 
#>                                                                                                                           bsdb:31976177/1/1_Obesity:obese-children_vs_normal-weight-children_UP 
#>                                                                                                                                                                                               5 
#>                                                                                                                           bsdb:31976177/2/1_Obesity:obese-children_vs_normal-weight-children_UP 
#>                                                                                                                                                                                              44 
#>                                                                                                                           bsdb:32143621/1/1_Obesity:Obese-children-(O)_vs_Normal-weight-(NW)_UP 
#>                                                                                                                                                                                               3 
#>                                                                                                        bsdb:32143621/2/1_Obesity:Obesity-with-metabolic-syndrome-(OMS)_vs_Normal-weight-(NW)_UP 
#>                                                                                                                                                                                              13 
#>                                          bsdb:32380969/1/1_Antimicrobial-agent:infants-who-received-short-term-or-long-term-antibiotic-treatment---week-1_vs_control-(no-treatment)---week-1_UP 
#>                                                                                                                                                                                               1 
#>                                                     bsdb:32380969/2/1_Antimicrobial-agent:infants-received-long-term(>5days)-antibiotic-treatment---week-2_vs_control(no-treatment)---week-2_UP 
#>                                                                                                                                                                                               1 
#>                                                                                                                                    bsdb:32429742/1/1_Obesity:obese_vs_normal-weight-children_UP 
#>                                                                                                                                                                                               1 
cn_res <- getCriticalN(relevant.sigs, my.siglengths)
cn_res$critical_n
#> 95% 
#>   9 
# Compare to observed
frequencySigs(my.sigs.increased)
#>   1578    976   1239   1591  39778    543   1350   1485 201174  31979 
#>      5      5      4      4      4      4      3      3      3      3 
```
