library(tidyverse)
library(scico)
library(ggpubr)
library(scales)
library(svglite)

## First two (Parts 1-2) sections of code create figures related to BOLD sequences only.
  ## Evaluations consist of BOLD data either untrimmed ("Full") or trimmed at particualr coordinates in sequence alignment ("ANML")
  ## BOLD data is also either dereplicated ('100') or clustered ('99 | 98 | 97'). See details below.

## Parts 3 and 4 evaluate ONLY the dereplicated, ANML-trimmed data, obtained from one of three sources:
  ## BOLD (as above, 'BOLD_ANML_100')
  ## NCBI-only-Bold ('ncbiOB'), obtained from GenBank requiring 'BARCODE' keyword... aka, cross referenced from BOLD dataset
  ## NCBI-no-Bold ('ncbiNB'), obtained from GenBank excluding 'BARCODE' keyword... aka, those that are not expected to be in BOLD dataset

################################################################################
## part 0a - import fitClassifier, crossValidate, evalSeqs, evalTaxa text files
################################################################################

## functions to import RESCRIPt outputs from `evaluate-seqs`, `evaluate-taxonomy`, `evaluate-cross-validate`, and `evaluate-fit-classifier`
  ## grouped by primer-trimmed (ANML) or not (Full), 
  ## and dereplication (100) or dereplication followed by clustering (97, 98, 99)
## datasets are divided into arthropod and chordate reference sequences


## import `evaluate-seqs` outputs, containing information from exporting all .qzv files into single text file per Taxa group (chordata, arthropoda)
evalSeqImporter <- function(urlpath, TaxaLabel){
  tmpSeqs <- read_delim(file=urlpath, col_names=FALSE, delim="\t")
  colnames(tmpSeqs) <- c("Metric", "Value", "Dataset")
  tmpSeqs <- tmpSeqs %>% 
    mutate(Dataset = gsub("boldDerep1", "boldFull", Dataset),
           TaxaLabel = TaxaLabel)
  tmpSeqs
}

evalSeq_chordateURL <- "https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper/data/evalSeqs/all_bold_ChordataOnly_EvalSeqs_Data.tsv"
evalSeq_chordate <- evalSeqImporter(evalSeq_chordateURL, "Chordata")
evalSeqs_arthropodURL <- "https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper/data/evalSeqs/all_bold_ArthOnly_EvalSeqs_Data.tsv"
evalSeq_arthropod <- evalSeqImporter(evalSeqs_arthropodURL, "Arthropoda")
evalSeqs <- rbind(evalSeq_chordate, evalSeq_arthropod)
rm(evalSeq_chordateURL, evalSeq_chordate, evalSeqs_arthropodURL, evalSeq_arthropod)


## import `evaluate-taxonomy` outputs, containing information from exporting all .qzv files into single text file per Taxa group (chordata, arthropoda)
evalTaxImporter <- function(urlpath, TaxaLabel){
  tmpTaxa <- read_delim(file=urlpath, col_names=FALSE, delim="\t")
  colnames(tmpTaxa) <- c("drop", "Level", "Unique Labels", "Taxonomic Entropy", "Terminal Labels", "Percent Terminal Labels","Lables at Depth",
                          "Fraction of Labels at Depth","Unclassified Labels", "Fraction of Unclassified Labels","Dataset")
  tmpTaxa %>% 
    select(-drop) %>% 
    mutate(Dataset = gsub("_Chordata", "", Dataset),
           Dataset = gsub("_Arthropoda", "", Dataset),
           TaxaLabel = TaxaLabel)
}

evalTax_chordateURL <- "https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper/data/evalTaxa/all_bold_Chordate_EvalTaxa_Data.tsv"
evalTax_chordate <- evalTaxImporter(evalTax_chordateURL, "Chordata")
evalTax_arthropodURL <- "https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper/data/evalTaxa/all_bold_Arth_EvalTaxa_data.tsv"
evalTax_arthropod <- evalTaxImporter(evalTax_arthropodURL, "Arthropoda")
evalTaxa <- rbind(evalTax_chordate, evalTax_arthropod)
rm(evalTax_chordate, evalTax_arthropod, evalTax_chordateURL, evalTax_arthropodURL)

## import `evaluate-fit-classifier` and `evaluate-cross-validate` outputs
evalOtherImporter <- function(urlpath, Taxalabel){
  tmpFitClass <- read_delim(file=urlpath, col_names=TRUE, delim="\t") %>% 
    mutate(Dataset = gsub("boldDerep1", "boldFull", Dataset),
           TaxaLabel = Taxalabel)
}

fitclass_chordateURL <- "https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper/data/fitClassifier/all_bold_Chordate_FitClassifier_data.tsv"
fitclass_chordate <- evalOtherImporter(fitclass_chordateURL, "Chordata")

fitclass_arthropodURL <- "https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper/data/fitClassifier/all_bold_Arth_fitclassifier_data.csv"
fitclass_arthropod <- read_csv(file=fitclass_arthropodURL) %>% mutate(TaxaLabel = "Arthropoda")
fitclass <- rbind(fitclass_chordate, fitclass_arthropod)
rm(fitclass_chordateURL, fitclass_chordate, fitclass_arthropodURL, fitclass_arthropod)

crossval_chordateURL <- "https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper/data/crossValidate/all_bold_Chordate_CrossValidate_data.tsv"
crossval_chordate <- evalOtherImporter(crossval_chordateURL, "Chordata")
crossval_arthropodURL <- "https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper/data/crossValidate/all_bold_Arth_crossvalidate_data.csv"
crossval_arthropod <- read_csv(file=crossval_arthropodURL) %>% mutate(TaxaLabel = "Arthropoda")
crossval <- rbind(crossval_chordate, crossval_arthropod)
rm(crossval_chordateURL, crossval_chordate, crossval_arthropodURL, crossval_arthropod)

################################################################################
## part 0b: get better colors for plots
################################################################################
ANMLpal4 <- scico(n = 4, begin = 0, end = 0.4, palette = "vik", direction = 1)
Fullpal4 <- scico(n = 4, begin = 0.6, end = 1, palette = "vik", direction = -1)
pal8 <- c(ANMLpal4, Fullpal4)


################################################################################
## part 1 - plot for Taxonomic information
## needs to be made separately, then stitched together
################################################################################

## parts 1a-c: Unique labels, Taxonomic Entropy, Terminal labels:
## reorder Dataset levels for legend:
evalTaxa$Dataset <- factor(evalTaxa$Dataset, levels = c(
  'boldFull_100', 'boldFull_99', 'boldFull_98', 'boldFull_97',
  'boldANML_100', 'boldANML_99', 'boldANML_98', 'boldANML_97'))

## plot 0 - use only for legend key:
p0 <- ggplot(evalTaxa, aes(x=Level, y=`Unique Labels`, color=Dataset)) +
  facet_wrap(~TaxaLabel, scales = "free_y") +
  geom_line(size=2) + theme_bw() + theme(panel.grid = element_blank(), 
                                         strip.text = element_text(size=12),
                                         legend.text = element_text(size=14),
                                         legend.title = element_text(size=16)) +
  scale_color_manual(values=pal8) + scale_y_continuous(labels=comma) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S")) +
  guides(color=guide_legend(nrow=4,byrow=FALSE))

## plot 1a:
p1a <- ggplot(evalTaxa, aes(x=Level, y=`Unique Labels`, color=Dataset)) +
  facet_wrap(~TaxaLabel, scales = "free_y") +
  geom_line() + theme_bw() + theme(panel.grid = element_blank(), strip.text = element_text(size=12), legend.position = "none") +
  scale_color_manual(values=pal8) + scale_y_continuous(labels=comma) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S")) +
  guides(color=guide_legend(nrow=4,byrow=FALSE))

## plot 1b:
p1b <- ggplot(evalTaxa, aes(x=Level, y=`Taxonomic Entropy`, color=Dataset)) +
  facet_wrap(~TaxaLabel) +
  geom_line() + theme_bw() + theme(panel.grid = element_blank(), strip.text = element_text(size=12), legend.position="none") +
  scale_color_manual(values=pal8) + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S"))

## plot 1c:
p1c <- ggplot(evalTaxa, aes(x=Level, y=`Terminal Labels`, color=Dataset)) +
  facet_wrap(~TaxaLabel, scales="free_y") +
  geom_line() + theme_bw() + theme(panel.grid = element_blank(), strip.text = element_text(size=12), legend.position="none") +
  scale_color_manual(values=pal8) + scale_y_continuous(labels=comma) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S"))

## plot 1d:
## specify order for legend
fitclass$Dataset <- factor(fitclass$Dataset, levels = c(
  'boldFull_100', 'boldFull_99', 'boldFull_98', 'boldFull_97',
  'boldANML_100', 'boldANML_99', 'boldANML_98', 'boldANML_97'))
## make plot
p1d <- ggplot(fitclass, aes(x=Level, y=`F-Measure`, color=Dataset)) +
  facet_wrap(~TaxaLabel) +
  geom_line() + theme_bw() + theme(panel.grid = element_blank(), strip.text = element_text(size=12), legend.position="none") +
  scale_color_manual(values=pal8) + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S"))

## plot 1e:
## specify order for legend

crossval$Dataset <- factor(crossval$Dataset, levels = c(
  'boldFull_100', 'boldFull_99', 'boldFull_98', 'boldFull_97',
  'boldANML_100', 'boldANML_99', 'boldANML_98', 'boldANML_97'))
## make plot
p1e <- ggplot(crossval, aes(x=Level, y=`F-Measure`, color=Dataset)) +
  facet_wrap(~TaxaLabel) +
  geom_line() + theme_bw() + theme(panel.grid = element_blank(), strip.text = element_text(size=12), legend.position="none") +
  scale_color_manual(values=pal8) + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S"))

## stitch together plots 1a-1e, grouping legend in bottom right of 3x2 grid:
commleg1 <- get_legend(p0)
ggarrange(p1a, p1b, p1c, p1d, commleg1, p1e,
          ncol=2, nrow=3, labels = c("A", "B", "C", "D", NA, "E"))

setwd("~/github/COIdatabases/rescript_paper/figures/")   ## change as needed
ggsave("boldCOI_TaxonomyEval_Figure.png", width=25, height=15, units="cm")
ggsave("boldCOI_TaxonomyEval_Figure.svg", width=25, height=15, units="cm")

## cleanup:
rm(p1a, p1b, p1c, p1d, p1e, commleg1, ANMLpal4, Fullpal4, p0)
   
################################################################################
## part 2 - plot for Sequence information
## also needs to be made separately, then stitched together, like with Plot1
################################################################################

## specify order for legend
evalSeqs$Dataset <- factor(evalSeqs$Dataset, levels = c(
  'boldFull_100', 'boldFull_99', 'boldFull_98', 'boldFull_97',
  'boldANML_100', 'boldANML_99', 'boldANML_98', 'boldANML_97'))

## plot 2a:
p2a <- ggplot(evalSeqs %>% filter(Metric=="N uniques"), 
       aes(x=Dataset, y=Value, fill=Dataset)) +
  facet_wrap(~TaxaLabel, scales = "free_y") +
  geom_col(position = position_dodge(width = 1.1)) +
  scale_y_continuous(labels=comma) +
  scale_fill_manual(values=pal8) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle=45, hjust=1),
        strip.text = element_text(size=12),
        legend.position = "top") +
  labs(x="", y="Unique Sequences") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

## plot 2b:
## subset dataset and set another group level for plot
kmerplotdat <- evalSeqs %>% filter(str_sub(Metric, start = -7) == "Entropy")
kmerplotdat$Metric <- gsub(" Entropy", "", kmerplotdat$Metric)
kmerplotdat$Metric <- factor(kmerplotdat$Metric, levels = c(
  "Sequence", "32mer", "24mer", "16mer", "8mer"))

p2b <- ggplot(kmerplotdat, aes(x=Metric, y=Value, color=Dataset, group=Dataset)) +
  facet_wrap(~TaxaLabel) +
  geom_line() +
  scale_color_manual(values=pal8) +
  theme_bw() + 
  theme(panel.grid = element_blank(), strip.text = element_text(size=12), legend.position = "none") +
  labs(x="", y="Entropy")


## stitch two together plots with 2 col, 1row grid:
ggarrange(p2a, p2b, nrow=2, ncol=1, 
          heights = c(1.8, 1),
          labels = c("A", "B"),
          align = "hv")

setwd("~/github/COIdatabases/rescript_paper/figures/")   ## change as needed
ggsave("boldCOI_SequenceEval_Figure.png", width=17, height=12, units="cm")
ggsave("boldCOI_SequenceEval_Figure.svg", width=17, height=12, units="cm")

## cleanup
rm(p2a, p2b)

################################################################################
## part 3 - import NCBI data, combine with ANML-only data from BOLD, and reformat for analysis
################################################################################

## Import NCBI taxa info, reformat labels, and combine with BOLD_ANML_100:
### first, import and reformat for evaluate-taxonomy data
ncbiTaxURL <- "https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper/data/evalTaxa/all_ncbi_EvalTaxa_data.tsv"
evalTaxa_ncbi <- read_delim(file=ncbiTaxURL, col_names = FALSE, delim="\t")
colnames(evalTaxa_ncbi) <- c("Level", "Unique Labels", "Taxonomic Entropy", "Terminal Labels", "Percent Terminal Labels","Lables at Depth",
                             "Fraction of Labels at Depth","Unclassified Labels", "Fraction of Unclassified Labels","Dataset")
evalTaxa_ncbi <- evalTaxa_ncbi %>% separate(Dataset, into=c("Dataset", "TaxaLabel"), sep="_")
evalTaxa_ncbi$TaxaLabel <- ifelse(evalTaxa_ncbi$TaxaLabel == "ArthOnly", "Arthropoda", "Chordata")
evalTaxa_boldSelect <- evalTaxa %>% 
  filter(Dataset == "boldANML_100") %>% 
  mutate(Dataset = gsub("_100", "", Dataset))
evalTaxa_NCBIandBOLD <- rbind(evalTaxa_boldSelect, evalTaxa_ncbi)


### second, import and reformat for fit classifier data
ncbi_fitClass_URL <- "https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper/data/fitClassifier/all_ncbi_fitclassifier_data.csv"
evalFC_ncbi <- read_csv(file=ncbi_fitClass_URL, col_names = TRUE)
evalFC_ncbi <- evalFC_ncbi %>% separate(Dataset, into=c("Dataset", "TaxaLabel"), sep="_")
evalFC_ncbi$TaxaLabel <- ifelse(evalFC_ncbi$TaxaLabel == "ArthOnly", "Arthropoda", "Chordata")
fitclass_boldSelect <- fitclass %>% filter(Dataset == "boldANML_100") %>% mutate(Dataset = gsub("_100", "", Dataset))
fitclass_NCBIandBOLD <- rbind(evalFC_ncbi, fitclass_boldSelect)

### third, import and reformat for cross validate data
ncbi_crossVal_URL <- 'https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper/data/crossValidate/all_ncbi_crossvalidate_data.csv'
evalCV_ncbi <- read_csv(file=ncbi_crossVal_URL, col_names = TRUE)
evalCV_ncbi <- evalCV_ncbi %>% separate(Dataset, into=c("Dataset", "TaxaLabel"), sep="_")
evalCV_ncbi$TaxaLabel <- ifelse(evalCV_ncbi$TaxaLabel == "ArthOnly", "Arthropoda", "Chordata")
crossval_boldSelect <- crossval %>% filter(Dataset == "boldANML_100") %>% mutate(Dataset = gsub("_100", "", Dataset))
crossval_NCBIandBOLD <- rbind(evalCV_ncbi, crossval_boldSelect)

### fourth, import and reformat evaluate-seq data
ncbi_evalSeqURL <- 'https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper/data/evalSeqs/all_ncbi_EvalSeqs_data.tsv'
evalSeqs_ncbi <- read_delim(file=ncbi_evalSeqURL, col_names=FALSE, delim="\t")
colnames(evalSeqs_ncbi) <- c("Metric", "Value", "Dataset")
evalSeqs_ncbi <- evalSeqs_ncbi %>% separate(Dataset, into=c("Dataset", "TaxaLabel"), sep="_")
evalSeqs_ncbi$TaxaLabel <- ifelse(evalSeqs_ncbi$TaxaLabel == "ArthOnly", "Arthropoda", "Chordata")  
evalSeqs_boldSelect <- evalSeqs %>% filter(Dataset == "boldANML_100") %>% mutate(Dataset = gsub("_100", "", Dataset))
evalSeqs_NCBIandBOLD <- rbind(evalSeqs_ncbi, evalSeqs_boldSelect)

## cleanup
rm(evalTaxa_ncbi, evalTaxa_boldSelect, evalTaxa, evalTaxImporter, ncbiTaxURL,
   evalFC_ncbi, fitclass_boldSelect, ncbi_fitClass_URL, evalOtherImporter, fitclass,
   crossval_boldSelect, evalCV_ncbi, ncbi_crossVal_URL, crossval,
   evalSeqs, evalSeqs_ncbi, evalSeqs_boldSelect, kmerplotdat, evalSeqImporter, ncbi_evalSeqURL)


################################################################################
## part4: plots
################################################################################

## first set of plots are for `evaluate-taxonomy`, `evaluate-fit-classifier`, and `evaluate-cross-validate` outputs
## second set of plots are for `evaluate-seqs` output

########################################
## taxonomy label plots
########################################

## custom color palette for NCBI, plus using original color code for 'boldANML_100':
## boldANML_100 == #601200
## ncbiOB == #a0af3a
## ncbiNB == #565d2f
## ncbiAll == #2a2c1b
#pal4 <- c('#601200', '#a0af3a', '#565d2f', '#2a2c1b')
pal4 <- c('#601200', '#2a2c1b', '#565d2f', '#a0af3a')

## specify order for legend
evalTaxa_NCBIandBOLD$Dataset <- factor(evalTaxa_NCBIandBOLD$Dataset, levels = c(
  "boldANML", "ncbiOB", "ncbiNB", "ncbiAll"))

## plot to get legend info
p3leg <- ggplot(evalTaxa_NCBIandBOLD, aes(x=Level, y=`Unique Labels`, color=Dataset)) +
  facet_wrap(~TaxaLabel, scales = "free_y") +
  geom_line(data=evalTaxa_NCBIandBOLD %>% filter(Dataset != "ncbiAll"), size=1.25, linetype="solid") + 
  geom_line(data=evalTaxa_NCBIandBOLD %>% filter(Dataset == "ncbiAll"), size=1.25, linetype="dashed") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), strip.text = element_text(size=12), 
        legend.text = element_text(size=14), legend.title = element_text(size=16)) +
  scale_color_manual(values=pal4) + scale_y_continuous(labels=comma) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S")) +
  guides(color=guide_legend(nrow=4,byrow=FALSE))


## plot 3a:
p3a <- ggplot(evalTaxa_NCBIandBOLD, aes(x=Level, y=`Unique Labels`, color=Dataset)) +
  facet_wrap(~TaxaLabel, scales = "free_y") +
  geom_line(data=evalTaxa_NCBIandBOLD %>% filter(Dataset != "ncbiAll"), size=0.75, linetype="solid") + 
  geom_line(data=evalTaxa_NCBIandBOLD %>% filter(Dataset == "ncbiAll"), size=0.5, linetype="dashed") + 
  theme_bw() + theme(panel.grid = element_blank(), strip.text = element_text(size=12), legend.position = "none") +
  scale_color_manual(values=pal4) + scale_y_continuous(labels=comma) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S")) +
  guides(color=guide_legend(nrow=4,byrow=FALSE))

## plot 3b:
p3b <- ggplot(evalTaxa_NCBIandBOLD, aes(x=Level, y=`Taxonomic Entropy`, color=Dataset)) +
  facet_wrap(~TaxaLabel) +
  geom_line(data=evalTaxa_NCBIandBOLD %>% filter(Dataset != "ncbiAll"), size=0.75, linetype="solid") + 
  geom_line(data=evalTaxa_NCBIandBOLD %>% filter(Dataset == "ncbiAll"), size=0.5, linetype="dashed") + 
  theme_bw() + theme(panel.grid = element_blank(), strip.text = element_text(size=12), legend.position="none") +
  scale_color_manual(values=pal4) + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S"))

## plot 3c:
p3c <- ggplot(evalTaxa_NCBIandBOLD, aes(x=Level, y=`Terminal Labels`, color=Dataset)) +
  facet_wrap(~TaxaLabel, scales="free_y") +
  geom_line(data=evalTaxa_NCBIandBOLD %>% filter(Dataset != "ncbiAll"), size=0.75, linetype="solid") + 
  geom_line(data=evalTaxa_NCBIandBOLD %>% filter(Dataset == "ncbiAll"), size=0.5, linetype="dashed") + 
  theme_bw() + theme(panel.grid = element_blank(), strip.text = element_text(size=12), legend.position="none") +
  scale_color_manual(values=pal4) + scale_y_continuous(labels=comma) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S"))

## plot 3d:
## specify order for legend
fitclass_NCBIandBOLD$Dataset <- factor(fitclass_NCBIandBOLD$Dataset, levels = c(
  "boldANML", "ncbiOB", "ncbiNB", "ncbiAll"))
## make plot
p3d <- ggplot(fitclass_NCBIandBOLD, aes(x=Level, y=`F-Measure`, color=Dataset)) +
  facet_wrap(~TaxaLabel) +
  geom_line(data=fitclass_NCBIandBOLD %>% filter(Dataset != "ncbiAll"), size=0.75, linetype="solid") + 
  geom_line(data=fitclass_NCBIandBOLD %>% filter(Dataset == "ncbiAll"), size=0.5, linetype="dashed") + 
  theme_bw() + theme(panel.grid = element_blank(), strip.text = element_text(size=12), legend.position="none") +
  scale_color_manual(values=pal4) + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S"))

## plot 3e:
## specify order for legend
crossval_NCBIandBOLD$Dataset <- factor(crossval_NCBIandBOLD$Dataset, levels = c(
  "boldANML", "ncbiOB", "ncbiNB", "ncbiAll"))
## make plot
p3e <- ggplot(crossval_NCBIandBOLD, aes(x=Level, y=`F-Measure`, color=Dataset)) +
  facet_wrap(~TaxaLabel) +
  geom_line(data=crossval_NCBIandBOLD %>% filter(Dataset != "ncbiAll"), size=0.75, linetype="solid") + 
  geom_line(data=crossval_NCBIandBOLD %>% filter(Dataset == "ncbiAll"), size=0.5, linetype="dashed") + 
  theme_bw() + theme(panel.grid = element_blank(), strip.text = element_text(size=12), legend.position="none") +
  scale_color_manual(values=pal4) + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S"))

## stitch together plots 1a-1e, grouping legend in bottom right of 3x2 grid:
commleg3 <- get_legend(p3leg)
ggarrange(p3a, p3b, p3c, p3d, commleg3, p3e,
          ncol=2, nrow=3, labels = c("A", "B", "C", "D", NA, "E"))

setwd("~/github/COIdatabases/rescript_paper/figures/")   ## change as needed
### need to modify legend to make 'ncbiAll' figure dashed instead of solid linetype
ggsave("ncbiANDboldCOI_TaxonomyEval_Figure.png", width=25, height=15, units="cm")
ggsave("ncbiANDboldCOI_TaxonomyEval_Figure.svg", width=25, height=15, units="cm")

## cleanup:
rm(p3a, p3b, p3c, p3d, p3e, p3leg, commleg3)


########################################
## sequence label plots
########################################

## plot 4a:
p4a <- ggplot(evalSeqs_NCBIandBOLD %>% filter(Metric=="N uniques"), 
              aes(x=Dataset, y=Value, fill=Dataset)) +
  facet_wrap(~TaxaLabel, scales = "free_y") +
  geom_col(position = position_dodge(width = 1.1)) +
  scale_y_continuous(labels=comma) +
  scale_fill_manual(values=pal4) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle=45, hjust=1),
        strip.text = element_text(size=12),
        legend.position = "top") +
  labs(x="", y="Unique Sequences") +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

## plot 4b:
## subset dataset and set another group level for plot
kmerplotdat_NCBIandBOLD <- evalSeqs_NCBIandBOLD %>% filter(str_sub(Metric, start = -7) == "Entropy")
kmerplotdat_NCBIandBOLD$Metric <- gsub(" Entropy", "", kmerplotdat_NCBIandBOLD$Metric)
kmerplotdat_NCBIandBOLD$Metric <- factor(kmerplotdat_NCBIandBOLD$Metric, levels = c(
  "Sequence", "32mer", "24mer", "16mer", "8mer"))

p4b <- ggplot(kmerplotdat_NCBIandBOLD, aes(x=Metric, y=Value, color=Dataset, group=Dataset)) +
  facet_wrap(~TaxaLabel) +
  geom_line(data=kmerplotdat_NCBIandBOLD %>% filter(Dataset != "ncbiAll"), size=0.75, linetype="solid") + 
  geom_line(data=kmerplotdat_NCBIandBOLD %>% filter(Dataset == "ncbiAll"), size=0.5, linetype="dashed") + 
  scale_color_manual(values=pal4) +
  theme_bw() + 
  theme(panel.grid = element_blank(), strip.text = element_text(size=12), legend.position = "none") +
  labs(x="", y="Entropy")


## stitch two together plots with 2 col, 1row grid:
ggarrange(p4a, p4b, nrow=2, ncol=1, 
          heights = c(1.8, 1),
          labels = c("A", "B"),
          align = "hv")

setwd("~/github/COIdatabases/rescript_paper/figures/")   ## change as needed
ggsave("ncbiANDboldCOI_SequenceEval_Figure.png", width=17, height=12, units="cm")
ggsave("ncbiANDboldCOI_SequenceEval_Figure.svg", width=17, height=12, units="cm")

## cleanup
rm(p4a, p4b)
