library(tidyverse)

################################################################################
## part 1 - import and reformat data
## NCBI data from `get-ncbi-data`
## BOLD data created as explained here: https://forum.qiime2.org/t/building-a-coi-database-from-bold-references/16129
################################################################################

########### (A) import ncbiOB data
NCBIob_df <- read_delim(file='https://github.com/devonorourke/COIdatabases/raw/master/rescript_paper/data/taxonomy_files/ncbiOB_taxonomy.tsv.gz',
                        delim="\t", col_names = TRUE)
## reformat retaining only Animalia taxa names, dropping redundant labels
## replace any spaces after first string in Class/Order/Family/Genus names
ncbiOB_tmp <- NCBIob_df %>%
  mutate(Taxon = gsub("k__Metazoa", "k__Animalia", Taxon)) %>%
  mutate(Taxon = gsub("; ", ";", Taxon)) %>%
  filter(grepl("k__Animalia",Taxon)) %>%
  select(Taxon) %>%
  distinct() %>%
  separate(Taxon, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  mutate(
    Class = gsub(" .*$", "", Class),
    Order = gsub(" .*$", "", Order),
    Family = gsub(" .*$", "", Family),
    Genus = gsub(" .*$", "", Genus)) %>%
  distinct() %>%
  mutate(RowNumber = row.names(.))
## remove instances where the same Class/Order/Family name is used (these are really not well described taxa)
ncbiOB_tmp_baddies <- ncbiOB_tmp %>% select(Class, Order, Family, RowNumber) %>%
  mutate(Class = gsub("c__", "", Class), Order = gsub("o__", "", Order), Family = gsub("f__", "", Family)) %>%
  mutate(TestCO = Class == Order,
         TestCF = Class == Family,
         TestOF = Order == Family) %>%
  filter(TestCO == TRUE & TestCF == TRUE & TestOF == TRUE)
## reformat to match BOLD references
ncbiOB_taxa <- ncbiOB_tmp %>%
  filter(!RowNumber %in% ncbiOB_tmp_baddies$RowNumber) %>%
  mutate(Class = gsub("c__Actinopterygii", "c__Actinopteri", Class),
         Class = gsub("c__Enopla", "c__Hoplonemertea", Class),
         Class = gsub("c__Hexapoda", "c__Diplura", Class),
         Class = gsub("c__Acoelomorpha", "c__Turbellaria", Class),
         Dataset = "ncbiOB") %>%
  select(-RowNumber)
## create species-specific reference list for our species-level comparisons:
ncbiOB_species <- ncbiOB_taxa %>%
  filter(Phylum != "p__" & Class != "c__" & Order != "o__" & Family != "f__" & Genus != "g__" & Species != "s__") %>%
  mutate(Species = gsub(" .*$", "", Species)) %>%
  filter(!grepl("environmental", Species)) %>%
  filter(!grepl("sp\\..*", Species)) %>%
  filter(!grepl("cf\\..*", Species)) %>%
  filter(!grepl("gen\\..*", Species)) %>%
  filter(!grepl("aff\\..*", Species)) %>%
  filter(!grepl("BOLD", Species)) %>%
  filter(!grepl("complex", Species)) %>%
  filter(!grepl("Janzen", Species)) %>%
  filter(!grepl("BIOUG", Species)) %>%
  distinct()
## cleanup
rm(NCBIob_df, ncbiOB_tmp, ncbiOB_tmp_baddies)


########### (B) import ncbiNB data
NCBInb_df <- read_delim(file='https://github.com/devonorourke/COIdatabases/raw/master/rescript_paper/data/taxonomy_files/ncbiNB_taxonomy.tsv.gz',
                        delim="\t", col_names = TRUE)
## reformat retaining only Animalia taxa names, dropping redundant labels
ncbiNB_tmp <- NCBInb_df %>%
  mutate(Taxon = gsub("k__Metazoa", "k__Animalia", Taxon)) %>%
  mutate(Taxon = gsub("; ", ";", Taxon)) %>%
  filter(grepl("k__Animalia",Taxon)) %>%
  select(Taxon) %>%
  distinct() %>%
  separate(Taxon, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  mutate(
    Class = gsub(" .*$", "", Class),
    Order = gsub(" .*$", "", Order),
    Family = gsub(" .*$", "", Family),
    Genus = gsub(" .*$", "", Genus)) %>%
  distinct() %>%
  mutate(RowNumber = row.names(.))
## identify instances where the same Class/Order/Family name is used (these are really not well described taxa)
ncbiNB_tmp_baddies <- ncbiNB_tmp %>% select(Class, Order, Family, RowNumber) %>%
  mutate(Class = gsub("c__", "", Class), Order = gsub("o__", "", Order), Family = gsub("f__", "", Family)) %>%
  mutate(TestCO = Class == Order,
         TestCF = Class == Family,
         TestOF = Order == Family) %>%
  filter(TestCO == TRUE & TestCF == TRUE & TestOF == TRUE)
## reformat to match BOLD references
ncbiNB_taxa <- ncbiNB_tmp %>%
  filter(!RowNumber %in% ncbiNB_tmp_baddies$RowNumber) %>%
  mutate(Class = gsub("c__Actinopterygii", "c__Actinopteri", Class),
         Class = gsub("c__Enopla", "c__Hoplonemertea", Class),
         Class = gsub("c__Hexapoda", "c__Diplura", Class),
         Class = gsub("c__Acoelomorpha", "c__Turbellaria", Class),
         Dataset = "ncbiNB") %>%
  select(-RowNumber)
## create species-specific reference list for our species-level comparisons:
ncbiNB_species <- ncbiNB_taxa %>%
  filter(Phylum != "p__" & Class != "c__" & Order != "o__" & Family != "f__" & Genus != "g__" & Species != "s__") %>%
  mutate(Species = gsub(" .*$", "", Species)) %>%
  filter(!grepl("environmental", Species)) %>%
  filter(!grepl("sp\\..*", Species)) %>%
  filter(!grepl("cf\\..*", Species)) %>%
  filter(!grepl("gen\\..*", Species)) %>%
  filter(!grepl("aff\\..*", Species)) %>%
  filter(!grepl("BOLD", Species)) %>%
  filter(!grepl("complex", Species)) %>%
  filter(!grepl("Janzen", Species)) %>%
  filter(!grepl("BIOUG", Species)) %>%
  distinct()
## cleanup
rm(NCBInb_df, ncbiNB_tmp, ncbiNB_tmp_baddies)


########### (C) import BOLD data
BOLD_df <- read_delim(file='https://github.com/devonorourke/COIdatabases/raw/master/rescript_paper/data/taxonomy_files/bold_taxonomy.tsv.gz',
                        delim="\t", col_names = TRUE)
## reformat to match NCBI references, retaining only Animalia taxa names, dropping redundant labels
## separate into taxonomic levels, dropping any Class names with uncertain taxa status
bold_taxa <- BOLD_df %>%
  mutate(Taxon = gsub("tax=", "", Taxon)) %>%
  filter(grepl("k__Animalia",Taxon)) %>%
  filter(!grepl("incertae_sedis",Taxon), ignore.case = TRUE) %>%
  filter(!grepl("Incertae_sedis",Taxon), ignore.case = TRUE) %>%
  select(Taxon) %>%
  distinct() %>%
  separate(Taxon, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
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
         Genus = gsub(" .*$", "", Genus)) %>%
  distinct()

## cleanup
rm(BOLD_df)


## combine all three datasets:
all_taxa <- rbind(bold_taxa, ncbiNB_taxa, ncbiOB_taxa)

################################################################################
## part 2 - gather values for 3 and 2-way intersections and unique taxa labels
## apply this for each taxa level from Phylum to Genus
################################################################################


################################################################################
## first for Phylum level
################################################################################

### 1. get counts of all samples with Phylum taxa for BOLD
nsamples_universe_BOLD_P <- bold_taxa %>% filter(Phylum != "p__") %>% nrow()
list_universe_BOLD_P <- bold_taxa %>% filter(Phylum != "p__") %>% distinct(Phylum) %>% pull()
ndistinct_BOLD_P <- length(list_universe_BOLD_P)
### 2. get the list of all possible uniuqe Phylum taxa for BOLD
nsamples_universe_ncbiOB_P <- ncbiOB_taxa %>% filter(Phylum != "p__") %>% nrow()
list_universe_ncbiOB_P <- ncbiOB_taxa %>% filter(Phylum != "p__") %>% distinct(Phylum) %>% pull()
ndistinct_ncbiOB_P <- length(list_universe_ncbiOB_P)
### 3. get the number of distinct Phylum taxa names from that list
nsamples_universe_ncbiNB_P <- ncbiNB_taxa %>% filter(Phylum != "p__") %>% nrow()
list_universe_ncbiNB_P <- ncbiNB_taxa %>% filter(Phylum != "p__") %>% distinct(Phylum) %>% pull()
ndistinct_ncbiNB_P <- length(list_universe_ncbiNB_P)
## what is the list of all possible Phylum names across all datasets?
### 4a. get a list of all possible distinct Phyla taxa across all 3 datasets (aka. outer join of 3 sets)
### 4b. how many distinct labels are there?
list_universe_all_P <- unique(c(list_universe_BOLD_P, list_universe_ncbiOB_P, list_universe_ncbiNB_P))
length_universe_all_P <- all_taxa %>% filter(Phylum %in% list_universe_all_P) %>% nrow()
### 5a. what Phylum names are shared across all 3 datasets? (aka. inner join of all 3 sets)
### 5b. get counts of these shared Phylum names
### 5c. what fraction of samples do these label represent in this group?
list_shared_all_P <- intersect(intersect(list_universe_BOLD_P, list_universe_ncbiOB_P), list_universe_ncbiNB_P)
length_shared_all_P <- all_taxa %>% filter(Phylum %in% list_shared_all_P) %>% nrow()
fraction_shared_all_P <- length_shared_all_P/length_universe_all_P  ## fraction shared is relative to sum of all 3 groups
### 6a. What Phylum names are shared across a pair of datasets, but not in the other? (aka. 2-set intersection)
### 6b. Get those counts
### 6c. what fraction of samples do these labels represent in this group?
## sets shorthand: set1==bold; set2==ncbiOB; set3==ncbiNB
## set12: found in bold and ncbiOB, but not ncbiNB
exclusive_set12_P <- setdiff(intersect(list_universe_BOLD_P, list_universe_ncbiOB_P), list_universe_ncbiNB_P)
length_exlusive_set12_P <- all_taxa %>% filter(Phylum %in% exclusive_set12_P) %>% nrow()
nsamples_universe_set12_P <- sum(nsamples_universe_BOLD_P, nsamples_universe_ncbiOB_P)
length_universe_all_P
exclusive_set13_P <- setdiff(intersect(list_universe_BOLD_P, list_universe_ncbiNB_P), list_universe_ncbiOB_P)
length_exlusive_set13_P <- all_taxa %>% filter(Phylum %in% exclusive_set13_P) %>% nrow()
fraction_exclusive_set13_P <- length_exlusive_set13_P / sum(nsamples_universe_BOLD_P, nsamples_universe_ncbiNB_P)
exclusive_set23_P <- setdiff(intersect(list_universe_ncbiOB_P, list_universe_ncbiNB_P), list_universe_BOLD_P)
length_exlusive_set23_P <- all_taxa %>% filter(Phylum %in% exclusive_set23_P) %>% nrow()
fraction_exclusive_set23_P <- length_exlusive_set23_P / sum(nsamples_universe_ncbiOB_P, nsamples_universe_ncbiNB_P)
### 7. get the unique stuff:
uniq_bold_P <- setdiff(list_universe_BOLD_P, c(list_universe_ncbiOB_P, list_universe_ncbiNB_P))
nsamples_unique_bold_P <- all_taxa %>% filter(Phylum %in% uniq_bold_P) %>% nrow()
fraction_unique_bold_P <- nsamples_unique_bold_P / nsamples_universe_BOLD_P   ## fraction is relative to it's own set size, not sum of all groups
uniq_ncbiOB_P <- setdiff(list_universe_ncbiOB_P, c(list_universe_BOLD_P, list_universe_ncbiNB_P))
nsamples_unique_ncbiOB_P <- all_taxa %>% filter(Phylum %in% uniq_ncbiOB_P) %>% nrow()
fraction_unique_ncbiOB_P <- nsamples_unique_ncbiOB_P / nsamples_universe_ncbiOB_P
uniq_ncbiNB_P <- setdiff(list_universe_ncbiNB_P, c(list_universe_BOLD_P, list_universe_ncbiOB_P))
nsamples_unique_ncbiNB_P <- all_taxa %>% filter(Phylum %in% uniq_ncbiNB_P) %>% nrow()
fraction_unique_ncbiNB_P <- nsamples_unique_ncbiNB_P / nsamples_universe_ncbiNB_P
### 8. Translate these values into a data.frame:
dat_P <- data.frame(
  TaxLevel = "Phylum",
  Set = c("Universe", "Set123", "Set12", "Set23", "Set13", "Set1", "Set2", "Set3"),
  nLabels = c(length(list_universe_all_P),length(list_shared_all_P),length(exclusive_set12_P),length(exclusive_set23_P),length(exclusive_set13_P),length(uniq_bold_P),length(uniq_ncbiOB_P),length(uniq_ncbiNB_P)),pSamples = c(1,fraction_shared_all_P,fraction_exclusive_set12_P,fraction_exclusive_set23_P,fraction_exclusive_set13_P,fraction_unique_bold_P,fraction_unique_ncbiOB_P,fraction_unique_ncbiNB_P))

