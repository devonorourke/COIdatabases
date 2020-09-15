library(tidyverse)
library(scales)
library(svglite)

################################################################################
## part 1 - import and reformat data
## will create datasets appropriate for Phylum through Genus-level analyses
## species-level analyses require much more specific filtering. see supplementary scripts at end of document for these
## NCBI data from `get-ncbi-data`
## BOLD data created as explained here: https://forum.qiime2.org/t/building-a-coi-database-from-bold-references/16129
################################################################################

##### NCBI import function:
ncbi_data_import_function <- function(URLpath, DatasetName){
  ## reformat retaining only Animalia taxa names, dropping redundant labels, and also...
  ## replace any spaces after first string in Class/Order/Family/Genus names
    raw_df <- read_delim(file=URLpath, delim="\t", col_names=TRUE) %>% 
    mutate(Taxon = gsub("k__Metazoa", "k__Animalia", Taxon)) %>%
    mutate(Taxon = gsub("; ", ";", Taxon)) %>%
    filter(grepl("k__Animalia",Taxon)) %>%
    select(Taxon) %>%
    separate(Taxon, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
    mutate(
      Class = gsub(" .*$", "", Class),
      Order = gsub(" .*$", "", Order),
      Family = gsub(" .*$", "", Family),
      Genus = gsub(" .*$", "", Genus)) %>%
    mutate(RowNumber = row.names(.))
  ## remove instances where the same Class/Order/Family name is used (these are really not well described taxa)
  tmp_baddies <- raw_df %>% select(Class, Order, Family, RowNumber) %>%
    mutate(Class = gsub("c__", "", Class), Order = gsub("o__", "", Order), Family = gsub("f__", "", Family)) %>%
    mutate(TestCO = Class == Order,
           TestCF = Class == Family,
           TestOF = Order == Family) %>%
    filter(TestCO == TRUE & TestCF == TRUE & TestOF == TRUE)
  ## reformat to match BOLD references
  ## remove instances where any Phylum through Genus-rank label contains a "."
  raw_df %>%
    filter(!RowNumber %in% tmp_baddies$RowNumber) %>%
    mutate(DotFilter = paste0(Phylum, Class, Order, Family, Genus)) %>% 
    filter(!grepl("\\.",DotFilter)) %>%
    mutate(Class = gsub("c__Actinopterygii", "c__Actinopteri", Class),
           Class = gsub("c__Enopla", "c__Hoplonemertea", Class),
           Class = gsub("c__Hexapoda", "c__Diplura", Class),
           Class = gsub("c__Acoelomorpha", "c__Turbellaria", Class),
           Dataset = DatasetName) %>%
    select(-RowNumber, -DotFilter)
}


########### (A) import ncbiOB data
ncbiOB_url <- 'https://github.com/devonorourke/COIdatabases/raw/master/rescript_paper/data/taxonomy_files/ncbiOB_taxonomy.tsv.gz'
ncbiOB_taxa <- ncbi_data_import_function(ncbiOB_url, "ncbiOB")

########### (B) import ncbiNB data
ncbiNB_url <- 'https://github.com/devonorourke/COIdatabases/raw/master/rescript_paper/data/taxonomy_files/ncbiNB_taxonomy.tsv.gz'
ncbiNB_taxa <- ncbi_data_import_function(ncbiNB_url, "ncbiNB")

########### (C) import BOLD data
bold_taxa <- read_delim(file='https://github.com/devonorourke/COIdatabases/raw/master/rescript_paper/data/taxonomy_files/bold_taxonomy.tsv.gz',
                      delim="\t", col_names = TRUE) %>%
  mutate(Taxon = gsub("tax=", "", Taxon)) %>%
  filter(grepl("k__Animalia",Taxon)) %>%
  filter(!grepl("incertae_sedis",Taxon), ignore.case = TRUE) %>%
  filter(!grepl("Incertae_sedis",Taxon), ignore.case = TRUE) %>%
  select(Taxon) %>%
  separate(Taxon, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  mutate(DotFilter = paste0(Phylum, Class, Order, Family, Genus)) %>% 
  filter(!grepl("\\.",DotFilter)) %>%
  select(-DotFilter) %>% 
  mutate(Class = gsub("c__Actinopterygii", "c__Actinopteri", Class),
         Class = gsub("c__Reptilia", "c__Lepidosauria", Class),
         Class = gsub("c__Elasmobranchii", "c__Chondrichthyes", Class),
         Class = gsub("c__Copepoda", "c__Hexanauplia", Class),
         Class = gsub("c__Thecostraca", "c__Hexanauplia", Class),
         Class = gsub("c__Enopla", "c__Hoplonemertea", Class),
         Class = gsub("c__Cephalaspidomorphi", "c__Hyperoartia", Class),
         Class = gsub("c__Holocephali", "c__Chondrichthyes", Class),
         Order = gsub("o__Diplura_order", "o__Diplura", Order),
         Dataset = "bold",
         Class = gsub(" .*$", "", Class),
         Order = gsub(" .*$", "", Order),
         Family = gsub(" .*$", "", Family),
         Genus = gsub(" .*$", "", Genus))


## combine all three datasets:
all_taxa <- rbind(bold_taxa, ncbiNB_taxa, ncbiOB_taxa)

rm(ncbiNB_taxa, ncbiOB_taxa, bold_taxa, ncbiNB_url, ncbiOB_url, ncbi_data_import_function)

################################################################################
## part 2 - gather values for 3 and 2-way intersections and unique taxa labels
## apply this for each taxa level from Phylum to Genus
## Species-level filtering requires additional considerations
################################################################################

#######################################
## first for Phylum level
#######################################

## gather the total length of all possible taxa labels, per Dataset
universe_df_P <- all_taxa %>% 
  filter(Phylum != "p__") %>% 
  group_by(Dataset, Phylum) %>% 
  summarise(nLabels = n()) %>% 
  ungroup() %>% 
  group_by(Dataset) %>% 
  summarise(sumLabels = sum(nLabels))
  
## gather data to tabulate intersections
intersection_df_P <- all_taxa %>% 
  filter(Phylum != "p__") %>% 
  group_by(Dataset, Phylum) %>% 
  summarise(nLabels = n()) %>% 
  pivot_wider(names_from = Dataset, values_from = nLabels)

## identify the complete cases where all three datasets contain the same taxa labels
all3_intersection_P <- intersection_df_P %>% filter(complete.cases(.))
## how many of those common taxa labels are shared across all 3 datasets?
all3_intersection_unique_taxaLabels_P <- data.frame(Value = nrow(all3_intersection_P),
                                                  Group = "allShared",
                                                  Rank = "Phylum")
## compare the original "universe" of all possible distinct sequences to those containing those common taxa labels shared in all 3 datasets
all3_intersection_nLabels_P <- all3_intersection_P %>% 
  pivot_longer(-Phylum, names_to="Dataset", values_to="nLabels") %>% 
  group_by(Dataset) %>% 
  summarise(nLabelsSet = sum(nLabels))
all3_intersection_nLabels_P <- merge(all3_intersection_nLabels_P, universe_df_P, by="Dataset") %>% 
  mutate(Rank = "Phylum", Group = "allShared")

## calculate the same values for the three possible 2-way intersections:
## set1 == bold; set2==ncbiOB; set3==ncbiNB
## here we're comparing 
## identify the three 2 way intersections
set12_intersection_data_P <- intersection_df_P %>% filter(!complete.cases(.)) %>% filter(bold > 0 & ncbiOB > 0 & is.na(ncbiNB))
set13_intersection_data_P <- intersection_df_P %>% filter(!complete.cases(.)) %>% filter(bold > 0 & ncbiNB > 0 & is.na(ncbiOB))
set23_intersection_data_P <- intersection_df_P %>% filter(!complete.cases(.)) %>% filter(ncbiOB > 0 & ncbiNB > 0 & is.na(bold))
## calculate how many common taxa labels are shared in each set
set12_intersection_unique_taxaLabels_P <- data.frame(Value = nrow(set12_intersection_data_P), Group = "set12", Rank = "Phylum")
set13_intersection_unique_taxaLabels_P <- data.frame(Value = nrow(set13_intersection_data_P), Group = "set13", Rank = "Phylum")
set23_intersection_unique_taxaLabels_P <- data.frame(Value = nrow(set23_intersection_data_P), Group = "set23", Rank = "Phylum")
## compare how many of these distinct sequences are in these intersecting sets
phylum_intersection_function <- function(inputData, GroupName){
  tmp_intersection_nLabels <- inputData %>% 
    pivot_longer(-Phylum, names_to="Dataset", values_to="nLabels") %>%   
    group_by(Dataset) %>% 
    summarise(nLabelsSet = sum(nLabels))
  merge(tmp_intersection_nLabels, universe_df_P, by="Dataset") %>% 
    mutate(Rank = "Phylum", Group = GroupName) %>% 
    filter(!is.na(nLabelsSet))
}
set12_intersection_nLabels_P <- phylum_intersection_function(set12_intersection_data_P, "set12")
set13_intersection_nLabels_P <- phylum_intersection_function(set13_intersection_data_P, "set13")
set23_intersection_nLabels_P <- phylum_intersection_function(set23_intersection_data_P, "set23")
  

## calculate the same values for the unique labels present in each of the datasets:
## identify the unique taxa labels
bold_uniq_data_P <- intersection_df_P %>% filter(!complete.cases(.)) %>% filter(bold > 0 & is.na(ncbiOB) & is.na(ncbiNB))
ncbiOB_uniq_data_P <- intersection_df_P %>% filter(!complete.cases(.)) %>% filter(ncbiOB > 0 & is.na(bold) & is.na(ncbiNB))
ncbiNB_uniq_data_P <- intersection_df_P %>% filter(!complete.cases(.)) %>% filter(ncbiNB > 0 & is.na(ncbiOB) & is.na(bold))
## calculate how many common taxa labels are shared in each set
bold_uniq_taxaLabels_P <- data.frame(Value = nrow(bold_uniq_data_P), Group = "uniq_bold", Rank = "Phylum")
ncbiOB_uniq_taxaLabels_P <- data.frame(Value = nrow(ncbiOB_uniq_data_P), Group = "uniq_ncbiOB", Rank = "Phylum")
ncbiNB_uniq_taxaLabels_P <- data.frame(Value = nrow(ncbiNB_uniq_data_P), Group = "uniq_ncbiNB", Rank = "Phylum")
## compare how many of these distinct sequences are in these unique sets
phylum_uniq_function <- function(inputData, GroupName){
  tmp_uniq_nLabels <- inputData %>% 
    pivot_longer(-Phylum, names_to="Dataset", values_to="nLabels") %>%   
    group_by(Dataset) %>% 
    summarise(nLabelsSet = sum(nLabels))
  merge(tmp_uniq_nLabels, universe_df_P, by="Dataset") %>% 
    mutate(Rank = "Phylum", Group = GroupName) %>% 
    filter(!is.na(nLabelsSet))
}
bold_uniq_nLabels_P <- phylum_intersection_function(bold_uniq_data_P, "uniq_bold")
ncbiOB_uniq_nLabels_P <- phylum_intersection_function(ncbiOB_uniq_data_P, "uniq_ncbiOB")
ncbiNB_uniq_nLabels_P <- phylum_intersection_function(ncbiNB_uniq_data_P, "uniq_ncbiNB")

## create data.frame collecting all data for label counts:
Phylum_nLabels <- rbind(all3_intersection_nLabels_P,
      set12_intersection_nLabels_P, set13_intersection_nLabels_P, set23_intersection_nLabels_P, 
      bold_uniq_nLabels_P, ncbiOB_uniq_nLabels_P, ncbiNB_uniq_nLabels_P)
## create data.frame collecting all data for unique taxa counts:
Phylum_nTaxa <- rbind(all3_intersection_unique_taxaLabels_P, 
                      set12_intersection_unique_taxaLabels_P, 
                      set13_intersection_unique_taxaLabels_P, 
                      set23_intersection_unique_taxaLabels_P, 
                      bold_uniq_taxaLabels_P, ncbiOB_uniq_taxaLabels_P, ncbiNB_uniq_taxaLabels_P)


## cleanup:
rm(list=ls(pattern = "*_P"))

#######################################
## for Class level
#######################################

## gather the total length of all possible taxa labels, per Dataset
universe_df_C <- all_taxa %>% 
  filter(Phylum != "p__" & Class != "c__") %>% 
  group_by(Dataset, Phylum, Class) %>% 
  summarise(nLabels = n()) %>% 
  ungroup() %>% 
  group_by(Dataset) %>% 
  summarise(sumLabels = sum(nLabels))

## gather data to tabulate intersections
intersection_df_C <- all_taxa %>% 
  filter(Phylum != "p__" & Class != "c__") %>% 
  mutate(Taxa = paste(Phylum,Class,sep=";")) %>% 
  group_by(Dataset, Taxa) %>% 
  summarise(nLabels = n()) %>% 
  pivot_wider(names_from = Dataset, values_from = nLabels)

## identify the complete cases where all three datasets contain the same taxa labels
all3_intersection_C <- intersection_df_C %>% filter(complete.cases(.))
## how many of those common taxa labels are shared across all 3 datasets?
all3_intersection_unique_taxaLabels_C <- data.frame(Value = nrow(all3_intersection_C),
                                                    Group = "allShared",
                                                    Rank = "Class")
## compare the original "universe" of all possible distinct sequences to those containing those common taxa labels shared in all 3 datasets
all3_intersection_nLabels_C <- all3_intersection_C %>% 
  pivot_longer(-Taxa, names_to="Dataset", values_to="nLabels") %>% 
  group_by(Dataset) %>% 
  summarise(nLabelsSet = sum(nLabels))
all3_intersection_nLabels_C <- merge(all3_intersection_nLabels_C, universe_df_C, by="Dataset") %>% 
  mutate(Rank = "Class", Group = "allShared")

## calculate the same values for the three possible 2-way intersections:
## set1 == bold; set2==ncbiOB; set3==ncbiNB
## here we're comparing 
## identify the three 2 way intersections
set12_intersection_data_C <- intersection_df_C %>% filter(!complete.cases(.)) %>% filter(bold > 0 & ncbiOB > 0 & is.na(ncbiNB))
set13_intersection_data_C <- intersection_df_C %>% filter(!complete.cases(.)) %>% filter(bold > 0 & ncbiNB > 0 & is.na(ncbiOB))
set23_intersection_data_C <- intersection_df_C %>% filter(!complete.cases(.)) %>% filter(ncbiOB > 0 & ncbiNB > 0 & is.na(bold))
## calculate how many common taxa labels are shared in each set
set12_intersection_unique_taxaLabels_C <- data.frame(Value = nrow(set12_intersection_data_C), Group = "set12", Rank = "Class")
set13_intersection_unique_taxaLabels_C <- data.frame(Value = nrow(set13_intersection_data_C), Group = "set13", Rank = "Class")
set23_intersection_unique_taxaLabels_C <- data.frame(Value = nrow(set23_intersection_data_C), Group = "set23", Rank = "Class")
## compare how many of these distinct sequences are in these intersecting sets
intersection_function_C <- function(inputData, GroupName){
  tmp_intersection_nLabels <- inputData %>% 
    pivot_longer(-Taxa, names_to="Dataset", values_to="nLabels") %>%   
    group_by(Dataset) %>% 
    summarise(nLabelsSet = sum(nLabels))
  merge(tmp_intersection_nLabels, universe_df_C, by="Dataset") %>% 
    mutate(Rank = "Class", Group = GroupName) %>% 
    filter(!is.na(nLabelsSet))
}
set12_intersection_nLabels_C <- intersection_function_C(set12_intersection_data_C, "set12")
set13_intersection_nLabels_C <- intersection_function_C(set13_intersection_data_C, "set13")
set23_intersection_nLabels_C <- intersection_function_C(set23_intersection_data_C, "set23")


## calculate the same values for the unique labels present in each of the datasets:
## identify the unique taxa labels
bold_uniq_data_C <- intersection_df_C %>% filter(!complete.cases(.)) %>% filter(bold > 0 & is.na(ncbiOB) & is.na(ncbiNB))
ncbiOB_uniq_data_C <- intersection_df_C %>% filter(!complete.cases(.)) %>% filter(ncbiOB > 0 & is.na(bold) & is.na(ncbiNB))
ncbiNB_uniq_data_C <- intersection_df_C %>% filter(!complete.cases(.)) %>% filter(ncbiNB > 0 & is.na(ncbiOB) & is.na(bold))
## calculate how many common taxa labels are shared in each set
bold_uniq_taxaLabels_C <- data.frame(Value = nrow(bold_uniq_data_C), Group = "uniq_bold", Rank = "Class")
ncbiOB_uniq_taxaLabels_C <- data.frame(Value = nrow(ncbiOB_uniq_data_C), Group = "uniq_ncbiOB", Rank = "Class")
ncbiNB_uniq_taxaLabels_C <- data.frame(Value = nrow(ncbiNB_uniq_data_C), Group = "uniq_ncbiNB", Rank = "Class")
## compare how many of these distinct sequences are in these unique sets
uniq_function_C <- function(inputData, GroupName){
  tmp_uniq_nLabels <- inputData %>% 
    pivot_longer(-Taxa, names_to="Dataset", values_to="nLabels") %>%   
    group_by(Dataset) %>% 
    summarise(nLabelsSet = sum(nLabels))
  merge(tmp_uniq_nLabels, universe_df_C, by="Dataset") %>% 
    mutate(Rank = "Class", Group = GroupName) %>% 
    filter(!is.na(nLabelsSet))
}
bold_uniq_nLabels_C <- uniq_function_C(bold_uniq_data_C, "uniq_bold")
ncbiOB_uniq_nLabels_C <- uniq_function_C(ncbiOB_uniq_data_C, "uniq_ncbiOB")
ncbiNB_uniq_nLabels_C <- uniq_function_C(ncbiNB_uniq_data_C, "uniq_ncbiNB")

## create data.frame collecting all data for label counts:
Class_nLabels <- rbind(all3_intersection_nLabels_C,
                        set12_intersection_nLabels_C, set13_intersection_nLabels_C, set23_intersection_nLabels_C, 
                        bold_uniq_nLabels_C, ncbiOB_uniq_nLabels_C, ncbiNB_uniq_nLabels_C)
## create data.frame collecting all data for unique taxa counts:
Class_nTaxa <- rbind(all3_intersection_unique_taxaLabels_C, 
                      set12_intersection_unique_taxaLabels_C, 
                      set13_intersection_unique_taxaLabels_C, 
                      set23_intersection_unique_taxaLabels_C, 
                      bold_uniq_taxaLabels_C, ncbiOB_uniq_taxaLabels_C, ncbiNB_uniq_taxaLabels_C)


## cleanup:
rm(list=ls(pattern = "*_C"))

#######################################
## for Order level
#######################################

## gather the total length of all possible taxa labels, per Dataset
universe_df_O <- all_taxa %>% 
  filter(Phylum != "p__" & Class != "c__" & Order != "o__") %>% 
  group_by(Dataset, Phylum, Class, Order) %>% 
  summarise(nLabels = n()) %>% 
  ungroup() %>% 
  group_by(Dataset) %>% 
  summarise(sumLabels = sum(nLabels))

## gather data to tabulate intersections
intersection_df_O <- all_taxa %>% 
  filter(Phylum != "p__" & Class != "c__" & Order != "o__") %>% 
  mutate(Taxa = paste(Phylum,Class,Order,sep=";")) %>% 
  group_by(Dataset, Taxa) %>% 
  summarise(nLabels = n()) %>% 
  pivot_wider(names_from = Dataset, values_from = nLabels)

## identify the complete cases where all three datasets contain the same taxa labels
all3_intersection_O <- intersection_df_O %>% filter(complete.cases(.))
## how many of those common taxa labels are shared across all 3 datasets?
all3_intersection_unique_taxaLabels_O <- data.frame(Value = nrow(all3_intersection_O),
                                                    Group = "allShared",
                                                    Rank = "Order")
## compare the original "universe" of all possible distinct sequences to those containing those common taxa labels shared in all 3 datasets
all3_intersection_nLabels_O <- all3_intersection_O %>% 
  pivot_longer(-Taxa, names_to="Dataset", values_to="nLabels") %>% 
  group_by(Dataset) %>% 
  summarise(nLabelsSet = sum(nLabels))
all3_intersection_nLabels_O <- merge(all3_intersection_nLabels_O, universe_df_O, by="Dataset") %>% 
  mutate(Rank = "Order", Group = "allShared")

## calculate the same values for the three possible 2-way intersections:
## set1 == bold; set2==ncbiOB; set3==ncbiNB
## here we're comparing 
## identify the three 2 way intersections
set12_intersection_data_O <- intersection_df_O %>% filter(!complete.cases(.)) %>% filter(bold > 0 & ncbiOB > 0 & is.na(ncbiNB))
set13_intersection_data_O <- intersection_df_O %>% filter(!complete.cases(.)) %>% filter(bold > 0 & ncbiNB > 0 & is.na(ncbiOB))
set23_intersection_data_O <- intersection_df_O %>% filter(!complete.cases(.)) %>% filter(ncbiOB > 0 & ncbiNB > 0 & is.na(bold))
## calculate how many common taxa labels are shared in each set
set12_intersection_unique_taxaLabels_O <- data.frame(Value = nrow(set12_intersection_data_O), Group = "set12", Rank = "Order")
set13_intersection_unique_taxaLabels_O <- data.frame(Value = nrow(set13_intersection_data_O), Group = "set13", Rank = "Order")
set23_intersection_unique_taxaLabels_O <- data.frame(Value = nrow(set23_intersection_data_O), Group = "set23", Rank = "Order")
## compare how many of these distinct sequences are in these intersecting sets
intersection_function_O <- function(inputData, GroupName){
  tmp_intersection_nLabels <- inputData %>% 
    pivot_longer(-Taxa, names_to="Dataset", values_to="nLabels") %>%   
    group_by(Dataset) %>% 
    summarise(nLabelsSet = sum(nLabels))
  merge(tmp_intersection_nLabels, universe_df_O, by="Dataset") %>% 
    mutate(Rank = "Order", Group = GroupName) %>% 
    filter(!is.na(nLabelsSet))
}
set12_intersection_nLabels_O <- intersection_function_O(set12_intersection_data_O, "set12")
set13_intersection_nLabels_O <- intersection_function_O(set13_intersection_data_O, "set13")
set23_intersection_nLabels_O <- intersection_function_O(set23_intersection_data_O, "set23")


## calculate the same values for the unique labels present in each of the datasets:
## identify the unique taxa labels
bold_uniq_data_O <- intersection_df_O %>% filter(!complete.cases(.)) %>% filter(bold > 0 & is.na(ncbiOB) & is.na(ncbiNB))
ncbiOB_uniq_data_O <- intersection_df_O %>% filter(!complete.cases(.)) %>% filter(ncbiOB > 0 & is.na(bold) & is.na(ncbiNB))
ncbiNB_uniq_data_O <- intersection_df_O %>% filter(!complete.cases(.)) %>% filter(ncbiNB > 0 & is.na(ncbiOB) & is.na(bold))
## calculate how many common taxa labels are shared in each set
bold_uniq_taxaLabels_O <- data.frame(Value = nrow(bold_uniq_data_O), Group = "uniq_bold", Rank = "Order")
ncbiOB_uniq_taxaLabels_O <- data.frame(Value = nrow(ncbiOB_uniq_data_O), Group = "uniq_ncbiOB", Rank = "Order")
ncbiNB_uniq_taxaLabels_O <- data.frame(Value = nrow(ncbiNB_uniq_data_O), Group = "uniq_ncbiNB", Rank = "Order")
## compare how many of these distinct sequences are in these unique sets
uniq_function_O <- function(inputData, GroupName){
  tmp_uniq_nLabels <- inputData %>% 
    pivot_longer(-Taxa, names_to="Dataset", values_to="nLabels") %>%   
    group_by(Dataset) %>% 
    summarise(nLabelsSet = sum(nLabels))
  merge(tmp_uniq_nLabels, universe_df_O, by="Dataset") %>% 
    mutate(Rank = "Order", Group = GroupName) %>% 
    filter(!is.na(nLabelsSet))
}
bold_uniq_nLabels_O <- uniq_function_O(bold_uniq_data_O, "uniq_bold")
ncbiOB_uniq_nLabels_O <- uniq_function_O(ncbiOB_uniq_data_O, "uniq_ncbiOB")
ncbiNB_uniq_nLabels_O <- uniq_function_O(ncbiNB_uniq_data_O, "uniq_ncbiNB")

## create data.frame collecting all data for label counts:
Order_nLabels <- rbind(all3_intersection_nLabels_O,
                       set12_intersection_nLabels_O, set13_intersection_nLabels_O, set23_intersection_nLabels_O, 
                       bold_uniq_nLabels_O, ncbiOB_uniq_nLabels_O, ncbiNB_uniq_nLabels_O)
## create data.frame collecting all data for unique taxa counts:
Order_nTaxa <- rbind(all3_intersection_unique_taxaLabels_O, 
                     set12_intersection_unique_taxaLabels_O, 
                     set13_intersection_unique_taxaLabels_O, 
                     set23_intersection_unique_taxaLabels_O, 
                     bold_uniq_taxaLabels_O, ncbiOB_uniq_taxaLabels_O, ncbiNB_uniq_taxaLabels_O)


## cleanup:
rm(list=ls(pattern = "*_O"))


#######################################
## for Family level
#######################################

## gather the total length of all possible taxa labels, per Dataset
universe_df_F <- all_taxa %>% 
  filter(Phylum != "p__" & Class != "c__" & Order != "o__" & Family != "f__") %>% 
  group_by(Dataset, Phylum, Class, Order, Family) %>% 
  summarise(nLabels = n()) %>% 
  ungroup() %>% 
  group_by(Dataset) %>% 
  summarise(sumLabels = sum(nLabels))

## gather data to tabulate intersections
intersection_df_F <- all_taxa %>% 
  filter(Phylum != "p__" & Class != "c__" & Order != "o__" & Family != "f__") %>% 
  mutate(Taxa = paste(Phylum,Class,Order,Family,sep=";")) %>% 
  group_by(Dataset, Taxa) %>% 
  summarise(nLabels = n()) %>% 
  pivot_wider(names_from = Dataset, values_from = nLabels)

## identify the complete cases where all three datasets contain the same taxa labels
all3_intersection_F <- intersection_df_F %>% filter(complete.cases(.))
## how many of those common taxa labels are shared across all 3 datasets?
all3_intersection_unique_taxaLabels_F <- data.frame(Value = nrow(all3_intersection_F),
                                                    Group = "allShared",
                                                    Rank = "Family")
## compare the original "universe" of all possible distinct sequences to those containing those common taxa labels shared in all 3 datasets
all3_intersection_nLabels_F <- all3_intersection_F %>% 
  pivot_longer(-Taxa, names_to="Dataset", values_to="nLabels") %>% 
  group_by(Dataset) %>% 
  summarise(nLabelsSet = sum(nLabels))
all3_intersection_nLabels_F <- merge(all3_intersection_nLabels_F, universe_df_F, by="Dataset") %>% 
  mutate(Rank = "Family", Group = "allShared")

## calculate the same values for the three possible 2-way intersections:
## set1 == bold; set2==ncbiOB; set3==ncbiNB
## here we're comparing 
## identify the three 2 way intersections
set12_intersection_data_F <- intersection_df_F %>% filter(!complete.cases(.)) %>% filter(bold > 0 & ncbiOB > 0 & is.na(ncbiNB))
set13_intersection_data_F <- intersection_df_F %>% filter(!complete.cases(.)) %>% filter(bold > 0 & ncbiNB > 0 & is.na(ncbiOB))
set23_intersection_data_F <- intersection_df_F %>% filter(!complete.cases(.)) %>% filter(ncbiOB > 0 & ncbiNB > 0 & is.na(bold))
## calculate how many common taxa labels are shared in each set
set12_intersection_unique_taxaLabels_F <- data.frame(Value = nrow(set12_intersection_data_F), Group = "set12", Rank = "Family")
set13_intersection_unique_taxaLabels_F <- data.frame(Value = nrow(set13_intersection_data_F), Group = "set13", Rank = "Family")
set23_intersection_unique_taxaLabels_F <- data.frame(Value = nrow(set23_intersection_data_F), Group = "set23", Rank = "Family")
## compare how many of these distinct sequences are in these intersecting sets
intersection_function_F <- function(inputData, GroupName){
  tmp_intersection_nLabels <- inputData %>% 
    pivot_longer(-Taxa, names_to="Dataset", values_to="nLabels") %>%   
    group_by(Dataset) %>% 
    summarise(nLabelsSet = sum(nLabels))
  merge(tmp_intersection_nLabels, universe_df_F, by="Dataset") %>% 
    mutate(Rank = "Family", Group = GroupName) %>% 
    filter(!is.na(nLabelsSet))
}
set12_intersection_nLabels_F <- intersection_function_F(set12_intersection_data_F, "set12")
set13_intersection_nLabels_F <- intersection_function_F(set13_intersection_data_F, "set13")
set23_intersection_nLabels_F <- intersection_function_F(set23_intersection_data_F, "set23")


## calculate the same values for the unique labels present in each of the datasets:
## identify the unique taxa labels
bold_uniq_data_F <- intersection_df_F %>% filter(!complete.cases(.)) %>% filter(bold > 0 & is.na(ncbiOB) & is.na(ncbiNB))
ncbiOB_uniq_data_F <- intersection_df_F %>% filter(!complete.cases(.)) %>% filter(ncbiOB > 0 & is.na(bold) & is.na(ncbiNB))
ncbiNB_uniq_data_F <- intersection_df_F %>% filter(!complete.cases(.)) %>% filter(ncbiNB > 0 & is.na(ncbiOB) & is.na(bold))
## calculate how many common taxa labels are shared in each set
bold_uniq_taxaLabels_F <- data.frame(Value = nrow(bold_uniq_data_F), Group = "uniq_bold", Rank = "Family")
ncbiOB_uniq_taxaLabels_F <- data.frame(Value = nrow(ncbiOB_uniq_data_F), Group = "uniq_ncbiOB", Rank = "Family")
ncbiNB_uniq_taxaLabels_F <- data.frame(Value = nrow(ncbiNB_uniq_data_F), Group = "uniq_ncbiNB", Rank = "Family")
## compare how many of these distinct sequences are in these unique sets
uniq_function_F <- function(inputData, GroupName){
  tmp_uniq_nLabels <- inputData %>% 
    pivot_longer(-Taxa, names_to="Dataset", values_to="nLabels") %>%   
    group_by(Dataset) %>% 
    summarise(nLabelsSet = sum(nLabels))
  merge(tmp_uniq_nLabels, universe_df_F, by="Dataset") %>% 
    mutate(Rank = "Family", Group = GroupName) %>% 
    filter(!is.na(nLabelsSet))
}
bold_uniq_nLabels_F <- uniq_function_F(bold_uniq_data_F, "uniq_bold")
ncbiOB_uniq_nLabels_F <- uniq_function_F(ncbiOB_uniq_data_F, "uniq_ncbiOB")
ncbiNB_uniq_nLabels_F <- uniq_function_F(ncbiNB_uniq_data_F, "uniq_ncbiNB")

## create data.frame collecting all data for label counts:
Family_nLabels <- rbind(all3_intersection_nLabels_F,
                       set12_intersection_nLabels_F, set13_intersection_nLabels_F, set23_intersection_nLabels_F, 
                       bold_uniq_nLabels_F, ncbiOB_uniq_nLabels_F, ncbiNB_uniq_nLabels_F)
## create data.frame collecting all data for unique taxa counts:
Family_nTaxa <- rbind(all3_intersection_unique_taxaLabels_F, 
                     set12_intersection_unique_taxaLabels_F, 
                     set13_intersection_unique_taxaLabels_F, 
                     set23_intersection_unique_taxaLabels_F, 
                     bold_uniq_taxaLabels_F, ncbiOB_uniq_taxaLabels_F, ncbiNB_uniq_taxaLabels_F)


## cleanup:
rm(list=ls(pattern = "*_F"))


#######################################
## for Genus level
#######################################

## gather the total length of all possible taxa labels, per Dataset
universe_df_G <- all_taxa %>% 
  filter(Phylum != "p__" & Class != "c__" & Order != "o__" & Family != "f__" & Genus != "g__") %>% 
  group_by(Dataset, Phylum, Class, Order, Family, Genus) %>% 
  summarise(nLabels = n()) %>% 
  ungroup() %>% 
  group_by(Dataset) %>% 
  summarise(sumLabels = sum(nLabels))

## gather data to tabulate intersections
intersection_df_G <- all_taxa %>% 
  filter(Phylum != "p__" & Class != "c__" & Order != "o__" & Family != "f__" & Genus != "g__") %>% 
  mutate(Taxa = paste(Phylum,Class,Order,Family,Genus,sep=";")) %>% 
  group_by(Dataset, Taxa) %>% 
  summarise(nLabels = n()) %>% 
  pivot_wider(names_from = Dataset, values_from = nLabels)

## identify the complete cases where all three datasets contain the same taxa labels
all3_intersection_G <- intersection_df_G %>% filter(complete.cases(.))
## how many of those common taxa labels are shared across all 3 datasets?
all3_intersection_unique_taxaLabels_G <- data.frame(Value = nrow(all3_intersection_G),
                                                    Group = "allShared",
                                                    Rank = "Genus")
## compare the original "universe" of all possible distinct sequences to those containing those common taxa labels shared in all 3 datasets
all3_intersection_nLabels_G <- all3_intersection_G %>% 
  pivot_longer(-Taxa, names_to="Dataset", values_to="nLabels") %>% 
  group_by(Dataset) %>% 
  summarise(nLabelsSet = sum(nLabels))
all3_intersection_nLabels_G <- merge(all3_intersection_nLabels_G, universe_df_G, by="Dataset") %>% 
  mutate(Rank = "Genus", Group = "allShared")

## calculate the same values for the three possible 2-way intersections:
## set1 == bold; set2==ncbiOB; set3==ncbiNB
## here we're comparing 
## identify the three 2 way intersections
set12_intersection_data_G <- intersection_df_G %>% filter(!complete.cases(.)) %>% filter(bold > 0 & ncbiOB > 0 & is.na(ncbiNB))
set13_intersection_data_G <- intersection_df_G %>% filter(!complete.cases(.)) %>% filter(bold > 0 & ncbiNB > 0 & is.na(ncbiOB))
set23_intersection_data_G <- intersection_df_G %>% filter(!complete.cases(.)) %>% filter(ncbiOB > 0 & ncbiNB > 0 & is.na(bold))
## calculate how many common taxa labels are shared in each set
set12_intersection_unique_taxaLabels_G <- data.frame(Value = nrow(set12_intersection_data_G), Group = "set12", Rank = "Genus")
set13_intersection_unique_taxaLabels_G <- data.frame(Value = nrow(set13_intersection_data_G), Group = "set13", Rank = "Genus")
set23_intersection_unique_taxaLabels_G <- data.frame(Value = nrow(set23_intersection_data_G), Group = "set23", Rank = "Genus")
## compare how many of these distinct sequences are in these intersecting sets
intersection_function_G <- function(inputData, GroupName){
  tmp_intersection_nLabels <- inputData %>% 
    pivot_longer(-Taxa, names_to="Dataset", values_to="nLabels") %>%   
    group_by(Dataset) %>% 
    summarise(nLabelsSet = sum(nLabels))
  merge(tmp_intersection_nLabels, universe_df_G, by="Dataset") %>% 
    mutate(Rank = "Genus", Group = GroupName) %>% 
    filter(!is.na(nLabelsSet))
}
set12_intersection_nLabels_G <- intersection_function_G(set12_intersection_data_G, "set12")
set13_intersection_nLabels_G <- intersection_function_G(set13_intersection_data_G, "set13")
set23_intersection_nLabels_G <- intersection_function_G(set23_intersection_data_G, "set23")


## calculate the same values for the unique labels present in each of the datasets:
## identify the unique taxa labels
bold_uniq_data_G <- intersection_df_G %>% filter(!complete.cases(.)) %>% filter(bold > 0 & is.na(ncbiOB) & is.na(ncbiNB))
ncbiOB_uniq_data_G <- intersection_df_G %>% filter(!complete.cases(.)) %>% filter(ncbiOB > 0 & is.na(bold) & is.na(ncbiNB))
ncbiNB_uniq_data_G <- intersection_df_G %>% filter(!complete.cases(.)) %>% filter(ncbiNB > 0 & is.na(ncbiOB) & is.na(bold))
## calculate how many common taxa labels are shared in each set
bold_uniq_taxaLabels_G <- data.frame(Value = nrow(bold_uniq_data_G), Group = "uniq_bold", Rank = "Genus")
ncbiOB_uniq_taxaLabels_G <- data.frame(Value = nrow(ncbiOB_uniq_data_G), Group = "uniq_ncbiOB", Rank = "Genus")
ncbiNB_uniq_taxaLabels_G <- data.frame(Value = nrow(ncbiNB_uniq_data_G), Group = "uniq_ncbiNB", Rank = "Genus")
## compare how many of these distinct sequences are in these unique sets
uniq_function_G <- function(inputData, GroupName){
  tmp_uniq_nLabels <- inputData %>% 
    pivot_longer(-Taxa, names_to="Dataset", values_to="nLabels") %>%   
    group_by(Dataset) %>% 
    summarise(nLabelsSet = sum(nLabels))
  merge(tmp_uniq_nLabels, universe_df_G, by="Dataset") %>% 
    mutate(Rank = "Genus", Group = GroupName) %>% 
    filter(!is.na(nLabelsSet))
}
bold_uniq_nLabels_G <- uniq_function_G(bold_uniq_data_G, "uniq_bold")
ncbiOB_uniq_nLabels_G <- uniq_function_G(ncbiOB_uniq_data_G, "uniq_ncbiOB")
ncbiNB_uniq_nLabels_G <- uniq_function_G(ncbiNB_uniq_data_G, "uniq_ncbiNB")

## create data.frame collecting all data for label counts:
Genus_nLabels <- rbind(all3_intersection_nLabels_G,
                        set12_intersection_nLabels_G, set13_intersection_nLabels_G, set23_intersection_nLabels_G, 
                        bold_uniq_nLabels_G, ncbiOB_uniq_nLabels_G, ncbiNB_uniq_nLabels_G)
## create data.frame collecting all data for unique taxa counts:
Genus_nTaxa <- rbind(all3_intersection_unique_taxaLabels_G, 
                      set12_intersection_unique_taxaLabels_G, 
                      set13_intersection_unique_taxaLabels_G, 
                      set23_intersection_unique_taxaLabels_G, 
                      bold_uniq_taxaLabels_G, ncbiOB_uniq_taxaLabels_G, ncbiNB_uniq_taxaLabels_G)


## cleanup:
rm(list=ls(pattern = "*_G"))


#######################################
## for Species level
#######################################

##### Species-level datasets require a bunch more filtering:
## filter NCBI data first
search_patterns <- c("(?=\\.*[A-Z])", "(?=.*\\d)")
ncbi_species_data_function <- function(NCBIdfName){
  NCBIdfName %>%
    filter(Phylum != "p__" & Class != "c__" & Order != "o__" & Family != "f__" & Genus != "g__" & Species != "s__") %>%
    mutate(Species = gsub(" .*$", "", Species)) %>%
    mutate(Species = gsub("\\(", "", Species)) %>%
    mutate(Species = gsub("\\)", "", Species)) %>%
    mutate(Species = gsub("\\'", "", Species)) %>%
    filter(!grepl("\\..*", Species)) %>%
    filter(!grepl("-", Species)) %>%
    filter(!grepl("environmental", Species)) %>%
    filter(!grepl("BOLD", Species)) %>%
    filter(!grepl("_complex$", Species)) %>%
    filter(!grepl("Janzen", Species)) %>%
    filter(!grepl("BIOUG", Species)) %>%
    filter(!grepl(paste0("(?=.*",search_patterns,")", collapse=""), Species, perl=TRUE)) %>%
    mutate(Species = tolower(Species))
}

ncbi_species <- ncbi_species_data_function(all_taxa %>% filter(Dataset != "bold"))

## filter BOLD data next
bold_species <- all_taxa %>% 
  filter(Dataset == "bold") %>% 
  filter(Phylum != "p__" & Class != "c__" & Order != "o__" & Family != "f__" & Genus != "g__" & Species != "s__") %>%
  filter(!grepl("\\..*", Species)) %>%
  filter(!grepl("-", Species)) %>%
  filter(!grepl("[0-9]", Species)) %>% 
  filter(!grepl("_complex$", Species)) %>%
  mutate(Species = gsub("\\(", "", Species)) %>%
  mutate(Species = gsub("\\)", "", Species)) %>%
  filter(!grepl(paste0("(?=.*",search_patterns,")", collapse=""), Species, perl=TRUE)) %>% 
  separate(Species, into = c("deleter", "Species"), sep=" ", extra = "drop", fill="right") %>% 
  select(-deleter) %>% 
  mutate(Species = paste0("s__", Species)) %>% 
  mutate(Species = tolower(Species))

## recombine:
all_species_taxa <- rbind(bold_species, ncbi_species)
## cleanup
rm(all_taxa, bold_species, ncbi_species)

## gather the total length of all possible taxa labels, per Dataset
universe_df_S <- all_species_taxa %>% 
  group_by(Dataset, Phylum, Class, Order, Family, Genus, Species) %>% 
  summarise(nLabels = n()) %>% 
  ungroup() %>% 
  group_by(Dataset) %>% 
  summarise(sumLabels = sum(nLabels))

## gather data to tabulate intersections
intersection_df_S <- all_species_taxa %>% 
  mutate(Taxa = paste(Phylum,Class,Order,Family,Genus,Species,sep=";")) %>% 
  group_by(Dataset, Taxa) %>% 
  summarise(nLabels = n()) %>% 
  pivot_wider(names_from = Dataset, values_from = nLabels)

## identify the complete cases where all three datasets contain the same taxa labels
all3_intersection_S <- intersection_df_S %>% filter(complete.cases(.))
## how many of those common taxa labels are shared across all 3 datasets?
all3_intersection_unique_taxaLabels_S <- data.frame(Value = nrow(all3_intersection_S),
                                                    Group = "allShared",
                                                    Rank = "Species")
## compare the original "universe" of all possible distinct sequences to those containing those common taxa labels shared in all 3 datasets
all3_intersection_nLabels_S <- all3_intersection_S %>% 
  pivot_longer(-Taxa, names_to="Dataset", values_to="nLabels") %>% 
  group_by(Dataset) %>% 
  summarise(nLabelsSet = sum(nLabels))
all3_intersection_nLabels_S <- merge(all3_intersection_nLabels_S, universe_df_S, by="Dataset") %>% 
  mutate(Rank = "Species", Group = "allShared")

## calculate the same values for the three possible 2-way intersections:
## set1 == bold; set2==ncbiOB; set3==ncbiNB
## here we're comparing 
## identify the three 2 way intersections
set12_intersection_data_S <- intersection_df_S %>% filter(!complete.cases(.)) %>% filter(bold > 0 & ncbiOB > 0 & is.na(ncbiNB))
set13_intersection_data_S <- intersection_df_S %>% filter(!complete.cases(.)) %>% filter(bold > 0 & ncbiNB > 0 & is.na(ncbiOB))
set23_intersection_data_S <- intersection_df_S %>% filter(!complete.cases(.)) %>% filter(ncbiOB > 0 & ncbiNB > 0 & is.na(bold))
## calculate how many common taxa labels are shared in each set
set12_intersection_unique_taxaLabels_S <- data.frame(Value = nrow(set12_intersection_data_S), Group = "set12", Rank = "Species")
set13_intersection_unique_taxaLabels_S <- data.frame(Value = nrow(set13_intersection_data_S), Group = "set13", Rank = "Species")
set23_intersection_unique_taxaLabels_S <- data.frame(Value = nrow(set23_intersection_data_S), Group = "set23", Rank = "Species")
## compare how many of these distinct sequences are in these intersecting sets
intersection_function_S <- function(inputData, GroupName){
  tmp_intersection_nLabels <- inputData %>% 
    pivot_longer(-Taxa, names_to="Dataset", values_to="nLabels") %>%   
    group_by(Dataset) %>% 
    summarise(nLabelsSet = sum(nLabels))
  merge(tmp_intersection_nLabels, universe_df_S, by="Dataset") %>% 
    mutate(Rank = "Species", Group = GroupName) %>% 
    filter(!is.na(nLabelsSet))
}
set12_intersection_nLabels_S <- intersection_function_S(set12_intersection_data_S, "set12")
set13_intersection_nLabels_S <- intersection_function_S(set13_intersection_data_S, "set13")
set23_intersection_nLabels_S <- intersection_function_S(set23_intersection_data_S, "set23")


## calculate the same values for the unique labels present in each of the datasets:
## identify the unique taxa labels
bold_uniq_data_S <- intersection_df_S %>% filter(!complete.cases(.)) %>% filter(bold > 0 & is.na(ncbiOB) & is.na(ncbiNB))
ncbiOB_uniq_data_S <- intersection_df_S %>% filter(!complete.cases(.)) %>% filter(ncbiOB > 0 & is.na(bold) & is.na(ncbiNB))
ncbiNB_uniq_data_S <- intersection_df_S %>% filter(!complete.cases(.)) %>% filter(ncbiNB > 0 & is.na(ncbiOB) & is.na(bold))
## calculate how many common taxa labels are shared in each set
bold_uniq_taxaLabels_S <- data.frame(Value = nrow(bold_uniq_data_S), Group = "uniq_bold", Rank = "Species")
ncbiOB_uniq_taxaLabels_S <- data.frame(Value = nrow(ncbiOB_uniq_data_S), Group = "uniq_ncbiOB", Rank = "Species")
ncbiNB_uniq_taxaLabels_S <- data.frame(Value = nrow(ncbiNB_uniq_data_S), Group = "uniq_ncbiNB", Rank = "Species")
## compare how many of these distinct sequences are in these unique sets
uniq_function_S <- function(inputData, GroupName){
  tmp_uniq_nLabels <- inputData %>% 
    pivot_longer(-Taxa, names_to="Dataset", values_to="nLabels") %>%   
    group_by(Dataset) %>% 
    summarise(nLabelsSet = sum(nLabels))
  merge(tmp_uniq_nLabels, universe_df_S, by="Dataset") %>% 
    mutate(Rank = "Species", Group = GroupName) %>% 
    filter(!is.na(nLabelsSet))
}
bold_uniq_nLabels_S <- uniq_function_S(bold_uniq_data_S, "uniq_bold")
ncbiOB_uniq_nLabels_S <- uniq_function_S(ncbiOB_uniq_data_S, "uniq_ncbiOB")
ncbiNB_uniq_nLabels_S <- uniq_function_S(ncbiNB_uniq_data_S, "uniq_ncbiNB")

## create data.frame collecting all data for label counts:
Species_nLabels <- rbind(all3_intersection_nLabels_S,
                         set12_intersection_nLabels_S, set13_intersection_nLabels_S, set23_intersection_nLabels_S, 
                         bold_uniq_nLabels_S, ncbiOB_uniq_nLabels_S, ncbiNB_uniq_nLabels_S)
## create data.frame collecting all data for unique taxa counts:
Species_nTaxa <- rbind(all3_intersection_unique_taxaLabels_S, 
                       set12_intersection_unique_taxaLabels_S, 
                       set13_intersection_unique_taxaLabels_S, 
                       set23_intersection_unique_taxaLabels_S, 
                       bold_uniq_taxaLabels_S, ncbiOB_uniq_taxaLabels_S, ncbiNB_uniq_taxaLabels_S)


## cleanup:
rm(list=ls(pattern = "*_S"))
rm(all_species_taxa, search_patterns, ncbi_species_data_function)

################################################################################
## part 3 - gather all output from part2 and plot
################################################################################

## gather data for both plot types
plotdat_nLabels <- rbind(Phylum_nLabels, Class_nLabels, Order_nLabels, Family_nLabels, Genus_nLabels, Species_nLabels)
#rm(Phylum_nLabels, Class_nLabels, Order_nLabels, Family_nLabels, Genus_nLabels, Species_nLabels)
plotdat_nTaxa <- rbind(Phylum_nTaxa, Class_nTaxa, Order_nTaxa, Family_nTaxa, Genus_nTaxa, Species_nTaxa)
#rm(Phylum_nTaxa, Class_nTaxa, Order_nTaxa, Family_nTaxa, Genus_nTaxa, Species_nTaxa)

## reorder labels for both plots
plotdat_nLabels$Rank <- factor(plotdat_nLabels$Rank, levels=c("Phylum", "Class", "Order", "Family", "Genus", "Species"))
plotdat_nTaxa$Rank <- factor(plotdat_nTaxa$Rank, levels=c("Phylum", "Class", "Order", "Family", "Genus", "Species"))


## plot1 - show the number of labels present among distinct sequences for a given taxa-rank (Phylum - Species)
ggplot(plotdat_nLabels, aes(group=Dataset, fill=Dataset)) +
  geom_col(aes(x=Group, y=sumLabels, group=Dataset),
           position = position_dodge2(preserve = "single"), inherit.aes = FALSE,
           fill = "white", color="black") +
  geom_col(aes(x=Group, y=nLabelsSet), position = position_dodge2(preserve = "single"), color="black") +
  facet_grid(Rank ~ .) +
  scale_y_continuous(labels = comma) +
  labs(x="", y="Sequences with taxon label") +
  scale_fill_manual(values = c("gray15", "gray50", "gray85")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

setwd("~/github/COIdatabases/rescript_paper/figures")
ggsave("COIlabelComps_noLabelValues_barStyle.svg", width=25, height=15, units="cm")  



## plot2 - show the number of labels present among distinct sequences for a given taxa-rank (Phylum - Species)
plot2_df <- plotdat_nTaxa %>% 
  pivot_wider(names_from = "Group", values_from="Value")
write_csv(plot2_df, path="comp_taxLabels_values.csv")

## make a quick data.table to write for the first plot, for discussion's sake:
plot1_df <- plotdat_nLabels %>% 
  mutate(fSeqsWLabel = nLabelsSet / sumLabels) %>% 
  mutate(fSeqsWLabel = round(fSeqsWLabel, 3)) %>% 
  select(Dataset, fSeqsWLabel, Rank, Group) %>% 
  pivot_wider(names_from = "Dataset", values_from="fSeqsWLabel") %>% 
  arrange(Rank, Group)

write_csv(plot1_df, path="comp_seqCount_forLabels.csv")
