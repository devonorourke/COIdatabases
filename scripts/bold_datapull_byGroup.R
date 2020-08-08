

## gathering Fungi, Protist, and Animal COI records from Barcode of Life Database Systems (BOLD)
## BOLD database queried 1-July-2020 (happy Canada Day)

## Need to update this script to check if new taxa names are used, as I manually ...
## ... copied the non-Arthropod, non-Chordate Animal names directly from BOLD website...
## ... http://v4.boldsystems.org/index.php/TaxBrowser_Home

library(bold)
library(taxize)
library(dplyr)


################################################################################
#### ---------------------- 0) The functions used ------------------------- ####
################################################################################

## filter bold data function:
gatherBOLDdat_function <- function(theboldlist){
  do.call(rbind.data.frame, theboldlist) %>%
  filter(markercode == "COI-5P") %>%
  select(sequenceID, processid, bin_uri, genbank_accession, nucleotides, country, institution_storing,
         phylum_name, class_name, order_name, family_name, genus_name, species_name)
}


## filtering a dataframe to get just the metadata function:
gatherBOLDmetadat_function <- function(thedataframe){
  thedataframe %>% 
    select(sequenceID, processid, bin_uri, genbank_accession, country, institution_storing)
}


## generate fasta function:
makefasta_function <- function(thedataframe, kingdomname){
  x.taxon <- thedataframe %>% select(sequenceID, phylum_name, class_name, order_name, family_name, genus_name, species_name)
  x.taxon$kingdom_name <- paste0("k__",kingdomname)
  x.taxon$phylum_name <- x.taxon$phylum_name %>% replace(is.na(.), "") %>% sub("^", "p__", .)
  x.taxon$class_name <- x.taxon$class_name %>% replace(is.na(.), "") %>% sub("^", "c__", .)
  x.taxon$order_name <- x.taxon$order_name %>% replace(is.na(.), "") %>% sub("^", "o__", .)
  x.taxon$family_name <- x.taxon$family_name %>% replace(is.na(.), "") %>% sub("^", "f__", .)
  x.taxon$genus_name <- x.taxon$genus_name %>% replace(is.na(.), "") %>% sub("^", "g__", .)
  x.taxon$species_name <- x.taxon$species_name %>% replace(is.na(.), "") %>% sub("^", "s__", .)
  x.taxon$taxon <- paste(x.taxon$kingdom_name, x.taxon$phylum_name, x.taxon$class_name,
                         x.taxon$order_name, x.taxon$family_name, x.taxon$genus_name,
                         x.taxon$species_name, sep = ";")
  x.taxon <- x.taxon %>% select(sequenceID, taxon)
  x.fasta <- thedataframe %>% select(sequenceID, nucleotides)
  merge(x.taxon, x.fasta)
}

################################################################################
#### ----- 1) get all the animals that aren't Chordates or Arthropods ----- ####
################################################################################

otherAnmlNames <- c("Acanthocephala","Acoelomorpha","Annelida","Brachiopoda","Bryozoa","Chaetognatha","Cnidaria","Ctenophora","Cycliophora","Echinodermata","Entoprocta","Gastrotricha","Gnathostomulida","Hemichordata","Kinorhyncha","Mollusca","Nematoda","Nematomorpha","Nemertea","Onychophora","Phoronida","Placozoa","Platyhelminthes","Porifera","Priapulida","Rhombozoa","Rotifera","Sipuncula","Tardigrada","Xenacoelomorpha")
altAnml_list <- lapply(otherAnmlNames, bold_seqspec)
altAnml_df <- gatherBOLDdat_function(altAnml_list)
alt_anml_fasta <- makefasta_function(altAnml_df, "Animalia")
write.csv(alt_anml_fasta, file = "boldCustom.allNonArthChordAnml.seqNtaxa.csv", quote = FALSE, row.names = FALSE)   ## file used to create taxonomy and fasta files
altAnml_meta <- gatherBOLDmetadat_function(altAnml_df)
write.csv(altAnml_meta, file = "boldCustom.allNonArthChordAnml.meta.csv", quote = FALSE, row.names = FALSE)

rm(altAnml_list, altAnml_df, alt_anml_fasta, altAnml_meta)

################################################################################
#### ------------------- 2) get all the fungi records --------------------- ####
################################################################################

fungiNames <- c('Ascomycota', 'Basidiomycota', 'Chytridiomycota', 'Glomeromycota', 'Myxomycota', 'Zygomycota')
fungi_list <- lapply(fungiNames, bold_seqspec)
fungi_df <- gatherBOLDdat_function(fungi_list)
fungi_fasta <- makefasta_function(fungi_df, "Fungi")
write.csv(fungi_fasta, file = "boldCustom.fungi.seqNtaxa.csv", quote = FALSE, row.names = FALSE)
fungi_meta <- gatherBOLDmetadat_function(fungi_df)
write.csv(fungi_meta, file = "boldCustom.fungi.meta.csv", quote = FALSE, row.names = FALSE)

rm(fungiNames, fungi_list, fungi_df, fungi_fasta, fungi_meta)

################################################################################
#### ------------------ 3) get all the protist records -------------------- ####
################################################################################

protistNames <- c('Chlorarachniophyta', 'Ciliophora', 'Heterokontophyta', 'Pyrrophycophyta')
protist_list <- lapply(protistNames, bold_seqspec)
protist_df <- gatherBOLDdat_function(protist_list)
protist_fasta <- makefasta_function(protist_df, "Protozoa")
write.csv(protist_fasta, file = "boldCustom.protist.seqNtaxa.csv", quote = FALSE, row.names = FALSE)
protist_meta <- gatherBOLDmetadat_function(protist_df)
write.csv(protist_meta, file = "boldCustom.protist.meta.csv", quote = FALSE, row.names = FALSE)

rm(protistNames, protist_list, protist_df, protist_fasta, protist_meta)

################################################################################
#### ----------------- 4) get all the chordate records -------------------- ####
################################################################################

chordateNames <- downstream("Chordata", db = "bold", downto = "class")
  # 14 classes represented. See 'chordate_list$Chordata$name'
chordate_list <- lapply(chordateNames$Chordata$name, bold_seqspec)
chordate_df <- gatherBOLDdat_function(chordate_list)
chordate_fasta <- makefasta_function(chordate_df, "Animalia")
write.csv(chordate_fasta, file = "boldCustom.chordate.seqNtaxa.csv", quote = FALSE, row.names = FALSE)
chordate_meta <- gatherBOLDmetadat_function(chordate_df)
write.csv(chordate_meta, file = "boldCustom.chordate.meta.csv", quote = FALSE, row.names = FALSE)

rm(chordateNames, chordate_list, chordate_df, chordate_fasta, chordate_meta)

################################################################################
#### Getting Arthropod records is tricky because we have to split up data into multiple groups:
#### All arthropods that are NOT insects
#### All insects that are NOT among the four largest orders (Coleoptera, Diptera, Hymenoptera, and Lepidoptera)
  #### Each of these four large orders (Col, Dip, Hymn, Lep) are further split too!
#### CAUTION: may need to modify these names as BOLD database grows
#### Haven't found a finite size where split needs to happen; generally > 500k records fail...
################################################################################

################################################################################
#### ---- 5) get all Arthropod records that are not of Insect Class ------- ####
################################################################################

allArthropod_names <- downstream("Arthropoda", db = "bold", downto = "class")
otherArth_names <- allArthropod_names$Arthropoda %>% filter(name != "Insecta") %>% select(name)
otherArth_list <- lapply(otherArth_names, bold_seqspec)
otherArth_df <- gatherBOLDdat_function(otherArth_list)
otherArth_fasta <- makefasta_function(otherArth_df, "Animalia")
write.csv(otherArth_fasta, file = "boldCustom.otherArthropods.seqNtaxa.csv", quote = FALSE, row.names = FALSE)
otherArth_meta <- gatherBOLDmetadat_function(otherArth_df)
write.csv(otherArth_meta, file = "boldCustom.otherArthropods.meta.csv", quote = FALSE, row.names = FALSE)

rm(otherArth_df, otherArth_fasta, otherArth_meta)

################################################################################
#### --- 6) get all Insect records for non-Col/Dip/Hym/Lep orders ------ ####
################################################################################

excludeNames <- c("Coleoptera","Diptera","Hymenoptera","Lepidoptera")
allInsect_names <- downstream("Insecta", db = "bold", downto = "order")
otherInsects_names <- allInsect_names$Insecta %>% filter(!name %in% excludeNames) %>% select(name)
otherInsects_list <- lapply(otherInsects_names, bold_seqspec)
otherInsects_df <- gatherBOLDdat_function(otherInsects_list)
otherInsects_fasta <- makefasta_function(otherInsects_df, "Animalia")
write.csv(otherInsects_fasta, file = "boldCustom.otherInsects.seqNtaxa.csv", quote = FALSE, row.names = FALSE)
otherInsects_meta <- gatherBOLDmetadat_function(otherInsects_df)
write.csv(otherInsects_meta, file = "boldCustom.otherInsects.meta.csv", quote = FALSE, row.names = FALSE)

rm(otherInsectOrder_names, otherInsects_list, otherInsects_df, otherInsects_fasta, otherInsects_meta)


################################################################################
#### ---------------- 7) get just Coleopteran records --------------------- ####
#### split into: Carabidae,Chrysomelidae,Curculionidae,Staphylinidae,(remaining others)
################################################################################

Col_list <- downstream("Coleoptera", db = "bold", downto = "family")
Col_Carabidae_names <- Col_list$Coleoptera %>% filter(name=="Carabidae") %>% select(name)
Col_Carabidae_list <- lapply(Col_Carabidae_names, bold_seqspec)
Col_Carabidae_df <- gatherBOLDdat_function(Col_Carabidae_list)
Col_Chrysomelidae_names <- Col_list$Coleoptera %>% filter(name=="Chrysomelidae") %>% select(name)
Col_Chrysomelidae_list <- lapply(Col_Chrysomelidae_names, bold_seqspec)
Col_Chrysomelidae_df <- gatherBOLDdat_function(Col_Chrysomelidae_list)
Col_Curculionidae_names <- Col_list$Coleoptera %>% filter(name=="Curculionidae") %>% select(name)
Col_Curculionidae_list <- lapply(Col_Curculionidae_names, bold_seqspec)
Col_Curculionidae_df <- gatherBOLDdat_function(Col_Curculionidae_list)
Col_Staphylinidae_names <- Col_list$Coleoptera %>% filter(name=="Staphylinidae") %>% select(name)
Col_Staphylinidae_list <- lapply(Col_Staphylinidae_names, bold_seqspec)
Col_Staphylinidae_df <- gatherBOLDdat_function(Col_Staphylinidae_list)
excludeColNames <- c("Carabidae","Chrysomelidae","Curculionidae","Staphylinidae")
Col_allother_names <- Col_list$Coleoptera %>% filter(!name %in% excludeColNames) %>% select(name)
Col_allothers_list <- lapply(Col_allother_names, bold_seqspec)
Col_allothers_df <- gatherBOLDdat_function(Col_allothers_list)
Col_df <- rbind(Col_Carabidae_df, Col_Chrysomelidae_df, Col_Curculionidae_df, Col_Staphylinidae_df, Col_allothers_df)   ## join all Coleopteran data
rm(Col_Carabidae_df, Col_Chrysomelidae_df, Col_Curculionidae_df, Col_Staphylinidae_df, Col_allothers_df)
rm(Col_Carabidae_names, Col_Chrysomelidae_names, Col_Curculionidae_names, Col_Staphylinidae_names, Col_allother_names)
rm(Col_Carabidae_list, Col_Chrysomelidae_list, Col_Curculionidae_list, Col_Staphylinidae_list, Col_allothers_list)
Col_fasta <- makefasta_function(Col_df, "Animalia")
write.csv(Col_fasta, file = "boldCustom.onlyColeoptera.seqNtaxa.csv", quote = FALSE, row.names = FALSE)
Col_meta <- gatherBOLDmetadat_function(Col_df)
write.csv(Col_meta, file = "boldCustom.onlyColeoptera.meta.csv", quote = FALSE, row.names = FALSE)
rm(Col_df, Col_fasta, Col_meta)

################################################################################
#### ------------------ 8) get just Dipteran records ---------------------- ####
#### split into: Sciaridae,Cecidomyiidae,Chironomidae,(remaining others)
################################################################################

diptera_list <- downstream("Diptera", db = "bold", downto = "family")
Dip_Sciaridae_names <- diptera_list$Diptera %>% filter(name=="Sciaridae") %>% select(name)
Dip_Sciaridae_list <- lapply(Dip_Sciaridae_names, bold_seqspec)
Dip_Sciaridae_df <- gatherBOLDdat_function(Dip_Sciaridae_list)
Dip_Cecidomyiidae_names <- diptera_list$Diptera %>% filter(name=="Cecidomyiidae") %>% select(name)
Dip_Cecidomyiidae_list <- lapply(Dip_Cecidomyiidae_names, bold_seqspec)
Dip_Cecidomyiidae_df <- gatherBOLDdat_function(Dip_Cecidomyiidae_list)
Dip_Chironomidae_names <- diptera_list$Diptera %>% filter(name=="Chironomidae") %>% select(name)
Dip_Chironomidae_list <- lapply(Dip_Chironomidae_names, bold_seqspec)
Dip_Chironomidae_df <- gatherBOLDdat_function(Dip_Chironomidae_list)
rm(Dip_Sciaridae_list, Dip_Cecidomyiidae_list, Dip_Chironomidae_list)
rm(Dip_Sciaridae_names, Dip_Cecidomyiidae_names, Dip_Chironomidae_list)
excludeDipNames <- c("Sciaridae","Cecidomyiidae","Chironomidae")
diptera_allother_names <- diptera_list$Diptera %>% filter(!name %in% excludeDipNames) %>% select(name)
Dip_allothers_list <- lapply(diptera_allother_names, bold_seqspec)
Dip_allothers_df <- gatherBOLDdat_function(Dip_allothers_list)
rm(Dip_Sciaridae_list, Dip_Cecidomyiidae_list, Dip_Chironomidae_list, Dip_allothers_list)
rm(Dip_Cecidomyiidae_names, Dip_Chironomidae_names, Dip_Sciaridae_names, diptera_allother_names)
Dip_df <- rbind(Dip_Sciaridae_df, Dip_Cecidomyiidae_df, Dip_Chironomidae_df, Dip_allothers_df)
rm(Dip_Sciaridae_df, Dip_Cecidomyiidae_df, Dip_Chironomidae_df, Dip_allothers_df)
Dip_fasta <- makefasta_function(Dip_df, "Animalia")
write.csv(Dip_fasta, file = "boldCustom.onlyDiptera.seqNtaxa.csv", quote = FALSE, row.names = FALSE)
Dip_meta <- gatherBOLDmetadat_function(Dip_df)
write.csv(Dip_meta, file = "boldCustom.onlyDiptera.meta.csv", quote = FALSE, row.names = FALSE)
rm(Dip_df, Dip_fasta, Dip_meta)


################################################################################
#### ---------------- 9) get just Hymenopteran records -------------------- ####
#### split into: Braconidae,Formicidae,Ichneumonidae,Platygastridae,(remaining others)
################################################################################

Hym_list <- downstream("Hymenoptera", db = "bold", downto = "family")
Hym_Braconidae_names <- Hym_list$Hymenoptera %>% filter(name=="Braconidae") %>% select(name)
Hym_Braconidae_list <- lapply(Hym_Braconidae_names, bold_seqspec)
Hym_Braconidae_df <- gatherBOLDdat_function(Hym_Braconidae_list)
Hym_Formicidae_names <- Hym_list$Hymenoptera %>% filter(name=="Formicidae") %>% select(name)
Hym_Formicidae_list <- lapply(Hym_Formicidae_names, bold_seqspec)
Hym_Formicidae_df <- gatherBOLDdat_function(Hym_Formicidae_list)
Hym_Ichneumonidae_names <- Hym_list$Hymenoptera %>% filter(name=="Ichneumonidae") %>% select(name)
Hym_Ichneumonidae_list <- lapply(Hym_Ichneumonidae_names, bold_seqspec)
Hym_Ichneumonidae_df <- gatherBOLDdat_function(Hym_Ichneumonidae_list)
Hym_Platygastridae_names <- Hym_list$Hymenoptera %>% filter(name=="Platygastridae") %>% select(name)
Hym_Platygastridae_list <- lapply(Hym_Platygastridae_names, bold_seqspec)
Hym_Platygastridae_df <- gatherBOLDdat_function(Hym_Platygastridae_list)
excludeColNames <- c("Braconidae","Formicidae","Ichneumonidae","Platygastridae")
Hym_allother_names <- Hym_list$Hymenoptera %>% filter(!name %in% excludeColNames) %>% select(name)
Hym_allothers_list <- lapply(Hym_allother_names, bold_seqspec)
Hym_allothers_df <- gatherBOLDdat_function(Hym_allothers_list)
rm(Hym_Braconidae_list, Hym_Formicidae_list, Hym_Ichneumonidae_list, Hym_Platygastridae_list, Hym_allothers_list)
rm(Hym_Braconidae_names, Hym_Formicidae_names, Hym_Ichneumonidae_names, Hym_Platygastridae_names, Hym_allothers_names)
Hym_df <- rbind(Hym_Braconidae_df, Hym_Formicidae_df, Hym_Ichneumonidae_df, Hym_Platygastridae_df, Hym_allothers_df)  ## join all Hymenopteran data
rm(Hym_Braconidae_df, Hym_Formicidae_df, Hym_Ichneumonidae_df, Hym_Platygastridae_df, Hym_allothers_df)
Hym_fasta <- makefasta_function(Hym_df, "Animalia")
write.csv(Hym_fasta, file = "boldCustom.onlyHymenoptera.seqNtaxa.csv", quote = FALSE, row.names = FALSE)
Hym_meta <- gatherBOLDmetadat_function(Hym_df)
write.csv(Hym_meta, file = "boldCustom.onlyHymenoptera.meta.csv", quote = FALSE, row.names = FALSE)
rm(Hym_df, Hym_fasta, Hym_meta)


################################################################################
#### ---------------- 10) get just Lepidoopteran records -------------------- ####
#### split into: Noctuidae,Erebidae,Sphingidae,Geometridae,(remaining others)
################################################################################

Lep_list <- downstream("Lepidoptera", db = "bold", downto = "family")
Lep_Noctuidae_names <- Lep_list$Lepidoptera %>% filter(name=="Noctuidae") %>% select(name)
Lep_Noctuidae_list <- lapply(Lep_Noctuidae_names, bold_seqspec)
Lep_Noctuidae_df <- gatherBOLDdat_function(Lep_Noctuidae_list)
Lep_Erebidae_names <- Lep_list$Lepidoptera %>% filter(name=="Erebidae") %>% select(name)
Lep_Erebidae_list <- lapply(Lep_Erebidae_names, bold_seqspec)
Lep_Erebidae_df <- gatherBOLDdat_function(Lep_Erebidae_list)
Lep_Sphingidae_names <- Lep_list$Lepidoptera %>% filter(name=="Sphingidae") %>% select(name)
Lep_Sphingidae_list <- lapply(Lep_Sphingidae_names, bold_seqspec)
Lep_Sphingidae_df <- gatherBOLDdat_function(Lep_Sphingidae_list)
Lep_Geometridae_names <- Lep_list$Lepidoptera %>% filter(name=="Geometridae") %>% select(name)
Lep_Geometridae_list <- lapply(Lep_Geometridae_names, bold_seqspec)
Lep_Geometridae_df <- gatherBOLDdat_function(Lep_Geometridae_list)
excludeColNames <- c("Noctuidae","Erebidae","Sphingidae","Geometridae")
Lep_allother_names <- Lep_list$Lepidoptera %>% filter(!name %in% excludeColNames) %>% select(name)
Lep_allothers_list <- lapply(Lep_allother_names, bold_seqspec)
Lep_allothers_df <- gatherBOLDdat_function(Lep_allothers_list)
rm(Lep_Noctuidae_names, Lep_Erebidae_names, Lep_Sphingidae_names, Lep_Geometridae_names, Lep_allother_names)
rm(Lep_Noctuidae_list, Lep_Erebidae_list, Lep_Sphingidae_list, Lep_Geometridae_list, Lep_allother_list)
Lep_df <- rbind(Lep_Noctuidae_df, Lep_Erebidae_df, Lep_Sphingidae_df, Lep_Geometridae_df, Lep_allothers_df)   ## join all Lepidopteran data
rm(Lep_Noctuidae_df, Lep_Erebidae_df, Lep_Sphingidae_df, Lep_Geometridae_df, Lep_allothers_df)
Lep_fasta <- makefasta_function(Lep_df, "Animalia")
write.csv(Lep_fasta, file = "boldCustom.onlyLepidoptera.seqNtaxa.csv", quote = FALSE, row.names = FALSE)
Lep_meta <- gatherBOLDmetadat_function(Lep_df)
write.csv(Lep_meta, file = "boldCustom.onlyLepidoptera.meta.csv", quote = FALSE, row.names = FALSE)
rm(Lep_df, Lep_fasta, Lep_meta)

