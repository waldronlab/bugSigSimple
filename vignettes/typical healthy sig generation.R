library(curatedMetagenomicData)
library(dplyr)
library(mia)
healthy <- sampleMetadata %>%
  filter(disease == 'healthy') %>%
  select(where( ~ !all(is.na(.x)))) 

adult <- healthy %>%
  filter(age_category %in% (c("adult", "senior")))

prev <- function(df, bodysite = "skin", prevalencecutoff = 0.000001) {
  .returnSig <- function(obj, threshold = 0, prevalence=0.000001){
    fractionpassing.logical <- rowSums(obj > threshold) / ncol(obj)
    rows.passing <- fractionpassing.logical > prevalence
    return(fractionpassing.logical[rows.passing])
  }
  df_sub <- df %>%
    filter(body_site == bodysite)
  df_sub1 <-
    df_sub %>% returnSamples("relative_abundance", rownames = "NCBI")
  df_sub2 <-
    df_sub %>% returnSamples("relative_abundance", rownames = "short")
  
  df_sub.byranks <- mia::splitByRanks(df_sub1)
  df_subg1 <-
    .returnSig(assay(df_sub.byranks[["genus"]]), prevalence = prevalencecutoff)
  df_subs1 <-
    .returnSig(assay(df_sub.byranks[["species"]]), prevalence = prevalencecutoff)
  
  df_sub.byranks <- mia::splitByRanks(df_sub2)
  df_subg2 <-
    .returnSig(assay(df_sub.byranks[["genus"]]), prevalence = prevalencecutoff)
  df_subs2 <-
    .returnSig(assay(df_sub.byranks[["species"]]), prevalence = prevalencecutoff)
  
  df_subg2 <- tibble(name = names(df_subg2), value = df_subg2)
  df_subs2 <- tibble(name = names(df_subs2), value = df_subs2)
  
  df_sub_species <- data.frame(df_subs1, df_subs2) %>%
    select(-c("df_subs1")) %>%
    arrange(desc(value))
  
  df_sub_genus <- data.frame(df_subg1, df_subg2) %>%
    select(-c("df_subg1")) %>%
    arrange(desc(value))
  return(list(species = df_sub_species, genus = df_sub_genus))
}

#adult
##skin
skinlist<-prev(adult, bodysite = "skin", prevalencecutoff = 0.000001)
skin_genus<-data.frame(skinlist[2])
skin_species<-data.frame(skinlist[1])

##vagina
vaginalist<-prev(adult, bodysite = "vagina", prevalencecutoff = 0.000001)
vagina_genus<-data.frame(vaginalist[2])
vagina_species<-data.frame(vaginalist[1])

##oralcavity
oralcavitylist<-prev(adult, bodysite = "oralcavity", prevalencecutoff = 0.000001)
oralcavity_genus<-data.frame(oralcavitylist[2])
oralcavity_species<-data.frame(oralcavitylist[1])

##nasalcavity
nasalcavitylist<-prev(adult, bodysite = "nasalcavity", prevalencecutoff = 0.000001)
nasalcavity_genus<-data.frame(nasalcavitylist[2])
nasalcavity_species<-data.frame(nasalcavitylist[1])

##stool
stoollist<-prev(adult, bodysite = "stool", prevalencecutoff = 0.000001)
stool_genus<-data.frame(stoollist[2])
stool_species<-data.frame(stoollist[1])

##milk
milklist<-prev(adult, bodysite = "milk", prevalencecutoff = 0.000001)
milk_genus<-data.frame(milklist[2])
milk_species<-data.frame(milklist[1])

## make 2 matrix- 1 for all species, 1 for all genus, stratified by body site
###genus
skin_genus$NCBI <- rownames(skin_genus)
vagina_genus$NCBI <- rownames(vagina_genus)
oralcavity_genus$NCBI <- rownames(oralcavity_genus)
nasalcavity_genus$NCBI <- rownames(nasalcavity_genus)
milk_genus$NCBI <- rownames(milk_genus)
stool_genus$NCBI <-rownames(stool_genus)

skin_genus<- skin_genus %>% dplyr::rename(skin_genus_prevalence=genus.value)
vagina_genus<- vagina_genus %>% dplyr::rename(vagina_genus_prevalence=genus.value)
oralcavity_genus<- oralcavity_genus %>% dplyr::rename(oral_genus_prevalence=genus.value)
nasalcavity_genus<- nasalcavity_genus %>% dplyr::rename(nasal_genus_prevalence=genus.value)
milk_genus<- milk_genus %>% dplyr::rename(milk_genus_prevalence=genus.value)
stool_genus<- stool_genus %>% dplyr::rename(stool_genus_prevalence=genus.value)

matrix_genus<- plyr::join_all(list(skin_genus,vagina_genus,oralcavity_genus,nasalcavity_genus,milk_genus,stool_genus), 
                              by = 'NCBI', type = 'full') 
matrix_genus[is.na(matrix_genus)] <- 0  
write.csv(matrix_genus,"~matrix_genus_adult.csv", row.names = TRUE)

###species
skin_species$NCBI <- rownames(skin_species)
vagina_species$NCBI <- rownames(vagina_species)
oralcavity_species$NCBI <- rownames(oralcavity_species)
nasalcavity_species$NCBI <- rownames(nasalcavity_species)
milk_species$NCBI <- rownames(milk_species)
stool_species$NCBI <-rownames(stool_species)

skin_species<- skin_species %>% dplyr::rename(skin_species_prevalence=species.value)
vagina_species<- vagina_species %>% dplyr::rename(vagina_species_prevalence=species.value)
oralcavity_species<- oralcavity_species %>% dplyr::rename(oral_species_prevalence=species.value)
nasalcavity_species<- nasalcavity_species %>% dplyr::rename(nasal_species_prevalence=species.value)
milk_species<- milk_species %>% dplyr::rename(milk_species_prevalence=species.value)
stool_species<- stool_species %>% dplyr::rename(stool_species_prevalence=species.value)

matrix_species<- plyr::join_all(list(skin_species,vagina_species,oralcavity_species,nasalcavity_species,milk_species,stool_species), 
                                by = 'NCBI', type = 'full') 
matrix_species[is.na(matrix_species)] <- 0  
write.csv(matrix_species,"~matrix_species_adult.csv", row.names = TRUE)

#child
##oralcavity
oralcavitylist<-prev(child, bodysite = "oralcavity", prevalencecutoff = 0.000001)
oralcavity_genusc<-data.frame(oralcavitylist[2])
oralcavity_speciesc<-data.frame(oralcavitylist[1])

##nasalcavity
nasalcavitylist<-prev(child, bodysite = "nasalcavity", prevalencecutoff = 0.000001)
nasalcavity_genusc<-data.frame(nasalcavitylist[2])
nasalcavity_speciesc<-data.frame(nasalcavitylist[1])

##stool
stoollist<-prev(child, bodysite = "stool", prevalencecutoff = 0.000001)
stool_genusc<-data.frame(stoollist[2])
stool_speciesc<-data.frame(stoollist[1])

## make 2 matrix- 1 for all species, 1 for all genus, stratified by body site
###genus
oralcavity_genusc$NCBI <- rownames(oralcavity_genusc)
nasalcavity_genusc$NCBI <- rownames(nasalcavity_genusc)
stool_genusc$NCBI <-rownames(stool_genusc)

oralcavity_genusc<- oralcavity_genusc %>% dplyr::rename(oral_genusc_prevalence=genus.value)
nasalcavity_genusc<- nasalcavity_genusc %>% dplyr::rename(nasal_genusc_prevalence=genus.value)
stool_genusc<- stool_genusc %>% dplyr::rename(stool_genusc_prevalence=genus.value)

matrix_genusc<- plyr::join_all(list(oralcavity_genusc,nasalcavity_genusc,stool_genusc), 
                               by = 'NCBI', type = 'full') 
matrix_genusc[is.na(matrix_genusc)] <- 0  
write.csv(matrix_genusc,"~matrix_genus_child.csv", row.names = TRUE)

###species
oralcavity_speciesc$NCBI <- rownames(oralcavity_speciesc)
nasalcavity_speciesc$NCBI <- rownames(nasalcavity_speciesc)
stool_speciesc$NCBI <-rownames(stool_speciesc)

oralcavity_speciesc<- oralcavity_speciesc %>% dplyr::rename(oral_speciesc_prevalence=species.value)
nasalcavity_speciesc<- nasalcavity_speciesc %>% dplyr::rename(nasal_speciesc_prevalence=species.value)
stool_speciesc<- stool_speciesc %>% dplyr::rename(stool_speciesc_prevalence=species.value)

matrix_speciesc<- plyr::join_all(list(oralcavity_speciesc,nasalcavity_speciesc,stool_speciesc), 
                                 by = 'NCBI', type = 'full') 
matrix_speciesc[is.na(matrix_speciesc)] <- 0  
write.csv(matrix_speciesc,"~matrix_species_child.csv", row.names = TRUE)