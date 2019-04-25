devtools::install_github("waldronlab/bugSigSimple")
library(bugSigSimple)
dat <- readCurationSheet("C:/Users/stath/Downloads/Microbial signatures curation - signatures (7).tsv")
dim(dat)
dat[1:5,1:5]


# print the column names of mydata2
names(dat)
names(dat)[23] <- "contrast"




my.contrasts<-c("taxons associated with BMI", "obese vs. controls", "obese adolescent vs. control adolescent", "obese adult vs. control adult",
                "bacterial taxa correlated with pediatric BMI", "pediatric obese vs. controls", "obese adolescent vs. controls",
                "pediatric obese vs. controls", "obese + metabolic syndrome vs. controls", "befores surgery vs. conrols", 
                "overweight/obese vs. controls", "abdominal obesity vs. controls")

my.pubmed.ids<-c("21829158", "302574444", "28628112", "30867711", "27499582",
                 "30568265", "25764541", "23032991", "23526699", "30669548",
                 "29388394", "29988340", "29950689", "30054529", "26230509",
                 "27007700", "27450202", "29922272", " 26261039", "27228093",
                 "29576948")
colnames(dat) 
rownames(dat)

ind1 <- dat[,"contrast"] %in% my.contrasts

my.dat <- dat[ind1 ,]

ind2 <- my.dat[, "PMID"] %in%  my.pubmed.ids

my.dat2  <- my.dat[,ind2]

my.dat3 <- subset(my.dat2, condition=="obesity")

dim(my.dat3)
