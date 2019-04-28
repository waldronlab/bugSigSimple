###devtools::install_github("waldronlab/bugSigSimple")
library(bugSigSimple)
dat <- readCurationSheet("C:/Users/stath/Downloads/Microbial signatures curation - signatures (9).tsv")
dim(dat)
dat[1:5,1:5]


# print the column names of mydata2
#names(dat)
names(dat)[23] <- "contrast"


sort(unique(dat$contrast))
dat$contrast[dat$contrast %in% c("obese adolescent vs. controls " ) ] <- "obese adolescent vs. controls"
dat$contrast[dat$contrast %in% c("obese vs. controls "  ) ] <- "obese vs. controls"
dat$contrast[dat$contrast %in% c("before surgery vs. controls "  ) ] <- "before surgery vs. controls"





my.contrasts<-c("taxons associated with BMI", "obese vs. controls", "obese adolescent vs. control adolescent ", "obese adult vs. control adult ",
                "obese adolescent vs. control adult " , "obese adult vs. control adolescent "  , "bacterial taxa correlated with pediatric BMI", "pediatric obese vs. controls",  "obese adolescent vs. controls",
                "pediatric obese vs. controls" , "obese + metabolic syndrome vs. controls", "before surgery vs. controls", 
                "overweight/obese vs. controls" , "abdominal obesity vs. controls", "Obesity group vs. controls")



my.pubmed.ids<-c("21829158", "302574444", "28628112", "30867711", "27499582",
                 "30568265", "25764541", "23032991", "23526699", "30669548",
                 "29388394", "29988340", "29950689", "30054529", "26230509",
                 "27007700", "27450202", "29922272", "26261039", "27228093",
                 "29576948")
colnames(dat) 
rownames(dat)

ind1 <- dat[,"contrast"] %in% my.contrasts
sum(ind1)
mean(ind1)
my.dat <- dat[ind1 ,]

ind2 <- my.dat[, "PMID"] %in%  my.pubmed.ids

my.dat2  <- my.dat[,ind2]

#my.dat3 <- subset(my.dat2, condition=="obesity")
my.dat3 <- subsetByCondition(my.dat2, condition="obesity")




getMostFrequentTaxa(my.dat3, n=40)
getMostFrequentTaxa(my.dat3, sig.type="UP", n=10)
getMostFrequentTaxa(my.dat3, sig.type="DOWN", n=10)


?binom.test
#binom test for frequent taxa "UP
binom.test(7, 12)
binom.test(6, 7)
binom.test(5, 5)
binom.test(5, 7)
binom.test(4, 4)
binom.test(3, 3)
binom.test(3, 4)
binom.test(3,3)

#binom test for frequent taxa "Down"
binom.test(5, 12)
binom.test(7, 7)
binom.test(4, 5)
binom.test(4, 4)
binom.test(3, 3)
binom.test(3, 4)
binom.test(2, 2)



?chisq.test

chisq.test(x=c(11,8,4,4), y=c(7-3,7-4,5-3,3-2), simulate.p.value = TRUE)

M <- as.table(rbind(c(2,5,0,0), c(11,8,4,4)- c(2,5,0,0)  ))
dimnames(M) <- list(down = c("Yes", "No"),
                    continent = c("Europe","Asia", "US", "Mexico"))
chisq.test(M, simulate.p.value = TRUE)

dim(my.dat3)

M <- as.table(rbind(c(2,5,0,0), c(7,7,4,3)- c(2,5,0,0)  ))
dimnames(M) <- list(down = c("Yes", "No"),
                    continent = c("Europe","Asia", "US", "Mexico"))
chisq.test(M, simulate.p.value = TRUE)

M


