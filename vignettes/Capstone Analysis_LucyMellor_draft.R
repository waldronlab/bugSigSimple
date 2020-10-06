#####################################
#                                   #
# Capstone                          #
#                                   #
# Lucy Mellor                       #
#####################################

#Install Packages
install.packages("devtools")
devtools::install_github("waldronlab/bugSigSimple")

#Load Library
library (bugSigSimple)
dat_all <- readCurationSheet("Microbial signatures curation - signatures (1).tsv")
dim(dat_all)


#Get subset of all signatures I curated or that were curated by others related my condition(s)
dat_subset <- subsetByCurator(dat_all, curator =  c("Lucy Mellor", "Rimsha Azhar"))
table(dat_subset[,"condition"])
dat_subset[,"condition"] <- sub(" $", "", dat_subset[,"condition"])

#Subset of conditions of interest
ind <- dat_subset[,"condition"] %in% c("asthma","atopic asthma","atopic eczema", "eczema", "egg allergy", "food allergy", "milk allergic reaction")
dat_subset2 <- dat_subset[ind,]
dim(dat_subset2)
table(dat_subset2[,"condition"])

#Name contrast column
names(dat_subset2)
names(dat_subset2)[26] <- "contrast"

sort(unique(dat_subset2$contrast))

#Remove studies that don't fit my inclusion criteria from analysis (i.e. publication year, adult, no healthy control group, non-human)
pubmed_ids <- c("10202341","22153774", "23339708", "27474122", "27838347", "29167295","29489859", "29518419", "29885665", "30643289", "31201890", "31619666")
ind <- dat_subset2[,"PMID"] %in% pubmed_ids
dat_subset3 <- dat_subset2[!ind, ]

#List of all contrasts
table(dat_subset3[,"contrast"])

#Remove contrasts that don't fit inclusion criteria from analysis (i.e. no healthy control group)
contrasts_controls <- c("infants IgE mediated vs. infants non-IgE mediated","atopic dermatitis after treatment with TCS + Bleach vs. baseline", "atopic dermatitis after treatment with TCS vs. baseline","healthy sibling vs. healthy control","7-18 years vs. <7 years","infant with non-IgE-mediated cow's milk allergy treated with EHCF vs. infant with non-IgE-mediated cow's milk allergy at diagnosis","infant with non-IgE-mediated cow's milk allergy treated with EHCF+LGG vs. infant with non-IgE-mediated cow's milk allergy at diagnosis","infant with non-IgE-mediated cow's milk allergy treated with EHCF+LGG vs. infant with non-IgE-mediated cow's milk allergy treated with EHCF","infant with non-IgE-mediated cow's milk allergy tolerant to cow's milk protein vs. infant with non-IgE-mediated cow's milk allergy non-tolerant to cow's milk protein","infants with persistent atopic dermatitis vs. transient atopic dermatitis", "Mothers whose infants developed IgE associated eczema vs. mothers whose infants remained non-allergic")
ind <- dat_subset3[,"contrast"] %in% contrasts_controls
dat_subset4 <- dat_subset3[!ind, ]

#List of included contrasts
table(dat_subset4[,"contrast"])

#List of included PMIDs
table(dat_subset4[,"PMID"])

#Frequencies of taxa in up and down (Overall)
getMostFrequentTaxa(dat_subset4, n=70)
getMostFrequentTaxa(dat_subset4,sig.type= "UP")
getMostFrequentTaxa(dat_subset4,sig.type= "DOWN")


#binomial test in UP
binom.test(x=7, n=7) #g__Faecalibacterium #p-value = 0.01563
binom.test(x=6, n=6) #g__Enterobacter #p-value = 0.03125
binom.test(x=5, n=5) #p__Bacteroidetes #p-value = 0.0625
binom.test(x=5, n=5) #c__Bacteroidia #p-value = 0.0625
binom.test(x=5, n=5) #o__Bacteroidales #p-value = 0.0625
binom.test(x=5, n=5) #g__Clostridium #p-value = 0.0625
binom.test(x=5, n=5) #g__Ruminococcus #p-value = 0.0625
binom.test(x=5, n=5) #g__Trabulsiella #p-value = 0.0625
binom.test(x=4, n=4) #g__Prevotella #p-value = 0.125
binom.test(x=4, n=4) #g__Staphylococcus p-value = 0.125

#binomial test in Down
binom.test(x=1, n=1) #g__Proteus #p-value = 1
binom.test(x=1, n=1) #g__Davidiella #p-value = 1
binom.test(x=1, n=1) #g__Sterigmatomyces #p-value = 1
binom.test(x=1, n=1) #g__Cryptococcus #p-value = 1

#Frequencies of taxa in up and down (Asthma)
ind1 <- dat_subset4[,"condition"] %in% c("asthma","atopic asthma")
dat_asthma <- dat_subset4[ind1,]
dim(dat_asthma)
table(dat_asthma[,"condition"])

getMostFrequentTaxa(dat_asthma, n=20)
getMostFrequentTaxa(dat_asthma,sig.type= "UP")
getMostFrequentTaxa(dat_asthma,sig.type= "DOWN")

#binomial test in UP
binom.test(x=3, n=3) #g_Prevotella #p-value = 0.25
binom.test(x=3, n=3) #g__Veillonella #p-value = 0.25
binom.test(x=1, n=1) #g__Rothia #p-value = 1
binom.test(x=1, n=1) #p__Firmicutes #p-value = 1
binom.test(x=1, n=1) #g__Staphylococcus #p-value = 1
binom.test(x=1, n=1) #g__Enterococcus #p-value = 1
binom.test(x=1, n=1) #g__Lactobacillus #p-value = 1
binom.test(x=1, n=1) #f__Clostridiaceae #p-value = 1
binom.test(x=1, n=1) #g__Clostridium #p-value = 1
binom.test(x=1, n=1) #g__Eubacterium #p-value = 1

#binomial test in Down
binom.test(x=1, n=1) #g__Proteus #p-value = 1
binom.test(x=1, n=1) #g__Davidiella #p-value = 1
binom.test(x=1, n=1) #g__Sterigmatomyces #p-value = 1
binom.test(x=1, n=1) #g__Cryptococcus #p-value = 1

#Frequencies of taxa in up and down (Food allergies)
ind2 <- dat_subset4[,"condition"] %in% c("egg allergy", "food allergy", "milk allergic reaction")
dat_FA <- dat_subset4[ind2,]
dim(dat_FA)
table(dat_FA[,"condition"])

getMostFrequentTaxa(dat_FA, n=20)
getMostFrequentTaxa(dat_FA,sig.type= "UP")
getMostFrequentTaxa(dat_FA,sig.type= "DOWN")

#binomial test in UP
binom.test(x=5, n=5) #p__Bacteroidetes #p-value = 0.0625
binom.test(x=5, n=5) #c__Bacteroidia #p-value = 0.0625
binom.test(x=5, n=5) #o__Bacteroidales #p-value = 0.0625
binom.test(x=5, n=5) #g__Enterobacter #p-value = 0.0625  
binom.test(x=5, n=5) #g__Trabulsiella #p-value = 0.0625
binom.test(x=4, n=4) #g__Clostridium #p-value = 0.125
binom.test(x=4, n=4) #g__Faecalibacterium #p-value = 0.125
binom.test(x=4, n=4) #g__Salmonella #p-value = 0.125
binom.test(x=3, n=3) #g__Ruminococcus #p-value = 0.25
binom.test(x=3, n=3) #f__Enterobacteriaceae #p-value = 0.25

#binomial test in Down
# none

#Frequencies of taxa in up and down (Eczema)
ind3 <- dat_subset4[,"condition"] %in% c("atopic eczema", "eczema")
dat_AD <- dat_subset4[ind3,]
dim(dat_AD)
table(dat_AD[,"condition"])

getMostFrequentTaxa(dat_AD, n=20)
getMostFrequentTaxa(dat_AD,sig.type= "UP")
getMostFrequentTaxa(dat_AD,sig.type= "DOWN")

#binomial test in UP
binom.test(x=3, n=3) #g__Parabacteroides #p-value = 0.25
binom.test(x=3, n=3) #g__Faecalibacterium #p-value = 0.25
binom.test(x=3, n=3) #g__Erwinia #p-value = 0.25
binom.test(x=2, n=2) #g__Bifidobacterium #p-value = 0.5
binom.test(x=2, n=2) #f__Porphyromonadaceae #p-value = 0.5
binom.test(x=2, n=2) #c__Bacilli #p-value = 0.5
binom.test(x=2, n=2) #g__Staphylococcus #p-value = 0.5
binom.test(x=2, n=2) #g__Streptococcus #p-value = 0.5
binom.test(x=2, n=2) #f__Lachnospiraceae #p-value = 0.5
binom.test(x=2, n=2) #g__Anaerostipes #p-value = 0.5

#binomial test in Down
# none

#List of body sites for included studies
table(dat_subset4[,"body_site"])

#Frequencies of taxa in up and down (feces)
ind4 <- dat_subset4[,"body_site"] %in% c("feces")
dat_feces <- dat_subset4[ind4,]
dim(dat_feces)
table(dat_feces[,"body_site"])

getMostFrequentTaxa(dat_feces, n=20)
getMostFrequentTaxa(dat_feces,sig.type= "UP")
getMostFrequentTaxa(dat_feces,sig.type= "DOWN")

#binomial test in UP
binom.test(x=7, n=7) #g__Faecalibacterium #p-value = 0.01563
binom.test(x=6, n=6) #g__Enterobacter #p-value = 0.03125
binom.test(x=5, n=5) #p__Bacteroidetes #p-value = 0.0625
binom.test(x=5, n=5) #c__Bacteroidia #p-value = 0.0625
binom.test(x=5, n=5) #o__Bacteroidales #p-value = 0.0625
binom.test(x=5, n=5) #g__Clostridium #p-value = 0.0625
binom.test(x=5, n=5) #g__Ruminococcus #p-value = 0.0625
binom.test(x=5, n=5) #g__Trabulsiella #p-value = 0.0625
binom.test(x=4, n=4) #f__Lachnospiraceae #p-value = 0.125
binom.test(x=4, n=4) #g__Roseburia #p-value = 0.125

#binomial test in Down
# none

#Frequencies of taxa in up and down (mouth/throat)
ind5 <- dat_subset4[,"body_site"] %in% c("hypopharynx", "nasopharynx", "oral gland", "throat")
dat_mouththroat <- dat_subset4[ind5,]
dim(dat_mouththroat)
table(dat_mouththroat[,"body_site"])

getMostFrequentTaxa(dat_mouththroat, n=20)
getMostFrequentTaxa(dat_mouththroat,sig.type= "UP")
getMostFrequentTaxa(dat_mouththroat,sig.type= "DOWN")

#binomial test in UP
binom.test(x=3, n=3) #g__Prevotella #p-value = 0.25
binom.test(x=1, n=1) #g__Staphylococcus #p-value = 1
binom.test(x=1, n=1) #g__Selenomonas #p-value = 1
binom.test(x=1, n=1) #g__Veillonella #p-value = 1
binom.test(x=1, n=1) #g__Parvimonas #p-value = 1
binom.test(x=1, n=1) #g__Fusobacterium #p-value = 1

#binomial test in Down
binom.test(x=1, n=1) #g__Proteus #p-value = 1
binom.test(x=1, n=1) #g__Davidiella #p-value = 1
binom.test(x=1, n=1) #g__Sterigmatomyces #p-value = 1
binom.test(x=1, n=1) #g__Cryptococcus #p-value = 1


#Frequencies of taxa in up and down (skin)
ind6 <- dat_subset4[,"body_site"] %in% c("skin of body", "skin of forearm")
dat_skin <- dat_subset4[ind6,]
dim(dat_skin)
table(dat_skin[,"body_site"])

getMostFrequentTaxa(dat_skin, n=20)
getMostFrequentTaxa(dat_skin,sig.type= "UP")
getMostFrequentTaxa(dat_skin,sig.type= "DOWN")

#binomial test in UP
binom.test(x=2, n=2) #g__Staphylococcus #p-value = 0.5
binom.test(x=1, n=1) #s__Corynebacterium_simulans #p-value = 1
binom.test(x=1, n=1) #g__Rhodococcus #p-value = 1
binom.test(x=1, n=1) #c__Bacilli #p-value = 1
binom.test(x=1, n=1) #s__Staphylococcus_aureus #p-value = 1
binom.test(x=1, n=1) #g__Alloiococcus #p-value = 1
binom.test(x=1, n=1) #s__Sphingomonas_yabuuchiae #p-value = 1
binom.test(x=1, n=1) #g__Tepidimonas #p-value = 1

#binomial test in Down
# none

binom.test(x=7, n=7) #g__Faecalibacterium #p-value = 0.01563
binom.test(x=6, n=6) #g__Enterobacter #p-value = 0.03125
binom.test(x=5, n=5) #p__Bacteroidetes #p-value = 0.0625
binom.test(x=5, n=5) #c__Bacteroidia #p-value = 0.0625
binom.test(x=5, n=5) #o__Bacteroidales #p-value = 0.0625
binom.test(x=5, n=5) #g__Clostridium #p-value = 0.0625
binom.test(x=5, n=5) #g__Ruminococcus #p-value = 0.0625
binom.test(x=5, n=5) #g__Trabulsiella #p-value = 0.0625
binom.test(x=4, n=4) #g__Prevotella #p-value = 0.125
binom.test(x=4, n=4) #g__Staphylococcus p-value = 0.125


#Frequencies of taxa in up and down (feces, asthma)
ind7 <- dat_asthma[,"body_site"] %in% c("feces")
dat_feces_asthma <- dat_asthma[ind7,]
dim(dat_feces_asthma)
table(dat_feces_asthma[,"body_site"])

getMostFrequentTaxa(dat_feces_asthma, n=20)
getMostFrequentTaxa(dat_feces_asthma,sig.type= "UP")
getMostFrequentTaxa(dat_feces_asthma,sig.type= "DOWN")

#Frequencies of taxa in up and down (feces, food allergies)
ind8 <- dat_FA[,"body_site"] %in% c("feces")
dat_feces_FA <- dat_FA[ind8,]
dim(dat_feces_FA)
table(dat_feces_FA[,"body_site"])

getMostFrequentTaxa(dat_feces_FA, n=20)
getMostFrequentTaxa(dat_feces_FA,sig.type= "UP")
getMostFrequentTaxa(dat_feces_FA,sig.type= "DOWN")

#Frequencies of taxa in up and down (feces, atopic dermatitis)
ind9 <- dat_AD[,"body_site"] %in% c("feces")
dat_feces_AD <- dat_FA[ind9,]
dim(dat_feces_AD)
table(dat_feces_AD[,"body_site"])

getMostFrequentTaxa(dat_feces_AD, n=20)
getMostFrequentTaxa(dat_feces_AD,sig.type= "UP")
getMostFrequentTaxa(dat_feces_AD,sig.type= "DOWN")
