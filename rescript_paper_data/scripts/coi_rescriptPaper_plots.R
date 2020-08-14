library(tidyverse)
library(scico)
library(ggpubr)

################################################################################
## part 0a - import fitClassifier, crossValidate, evalSeqs, evalTaxa text files
################################################################################

evalSeqs <- read_delim(file="https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper_data/all_Chordate_EvalSeqs_Data.tsv",
                       col_names = FALSE, delim="\t")
colnames(evalSeqs) <- c("Metric", "Value", "Dataset")
evalSeqs$Dataset <- gsub("boldDerep1", "boldFull", evalSeqs$Dataset)

evalTaxa <- read_delim(file="https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper_data/all_Chordate_EvalTaxa_Data.tsv",
                       col_names = FALSE, delim="\t") 
evalTaxa <- evalTaxa[,c(2:11)]
colnames(evalTaxa) <- c("Level", "Unique Labels", "Taxonomic Entropy",
                        "Terminal Labels", "Percent Terminal Labels",
                        "Lables at Depth", "Fraction of Labels at Depth",
                        "Unclassified Labels", "Fraction of Unclassified Labels",
                        "Dataset") 
evalTaxa$Dataset <- gsub("boldDerep1", "boldFull", evalTaxa$Dataset)
evalTaxa$Dataset <- gsub("_Chordata", "", evalTaxa$Dataset)

fitclass <- read_delim(file="https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper_data/all_Chordate_FitClassifier_Data.tsv",
                       col_names = TRUE, delim="\t") %>% 
  mutate(Dataset = gsub("boldDerep1", "boldFull", Dataset))


crossval <- read_delim(file="https://raw.githubusercontent.com/devonorourke/COIdatabases/master/rescript_paper_data/all_Chordate_CrossValidate_Data.tsv",
                       col_names = TRUE, delim="\t") %>% 
  mutate(Dataset = gsub("boldDerep1", "boldFull", Dataset))


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

## plot 1a:
p1a <- ggplot(evalTaxa, aes(x=Level, y=`Unique Labels`, color=Dataset)) +
  geom_line() + theme_bw() + theme(panel.grid = element_blank()) +
  scale_color_manual(values=pal8) + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S"))


## plot 1b:
p1b <- ggplot(evalTaxa, aes(x=Level, y=`Taxonomic Entropy`, color=Dataset)) +
  geom_line() + theme_bw() + theme(panel.grid = element_blank(), legend.position="none") +
  scale_color_manual(values=pal8) + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S"))


## plot 1c with alternate palette:
p1c <- ggplot(evalTaxa, aes(x=Level, y=`Terminal Labels`, color=Dataset)) +
  geom_line() + theme_bw() + theme(panel.grid = element_blank(), legend.position="none") +
  scale_color_manual(values=pal8) + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S"))

## plot 1d:
## specify order for legend
fitclass$Dataset <- factor(fitclass$Dataset, levels = c(
  'boldFull_100', 'boldFull_99', 'boldFull_98', 'boldFull_97',
  'boldANML_100', 'boldANML_99', 'boldANML_98', 'boldANML_97'))
## make plot
p1d <- ggplot(fitclass, aes(x=Level, y=`F-Measure`, color=Dataset)) +
  geom_line() + theme_bw() + theme(panel.grid = element_blank(), legend.position="none") +
  scale_color_manual(values=pal8) + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S"))

## plot 1e:
## specify order for legend
crossval$Dataset <- factor(crossval$Dataset, levels = c(
  'boldFull_100', 'boldFull_99', 'boldFull_98', 'boldFull_97',
  'boldANML_100', 'boldANML_99', 'boldANML_98', 'boldANML_97'))
## make plot
p1e <- ggplot(crossval, aes(x=Level, y=`F-Measure`, color=Dataset)) +
  geom_line() + theme_bw() + theme(panel.grid = element_blank(), legend.position="none") +
  scale_color_manual(values=pal8) + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("K", "P", "C", "O", "F", "G", "S"))

## stitch together plots 1a-1e, grouping legend in bottom right of 3x2 grid:
commleg1a <- get_legend(p1a)
mod_p1a <- p1a + theme(legend.position = "none")
ggarrange(mod_p1a, p1b, p1c, p1d, p1e, commleg1a,
          ncol=3, nrow=2, labels = c("A", "B", "C", "D", "E", NA))

setwd("~/github/COIdatabases/rescript_paper_data/")   ## change as needed
ggsave("boldCOI_chordateOnly_TaxonomyEval_Figure.png", width=25, height=15, units="cm")

## cleanup:
rm(p1a, mod_p1a, p1b, p1c, p1d, p1e, commleg1a,
   crossval, fitclass, evalTaxa, ANMLpal4, Fullpal4)
   
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
  geom_col(position = position_dodge(width = 1.1)) +
  scale_fill_manual(values=pal8) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="", y="Unique Sequences")


## plot 2b:
## subset dataset and set another group level for plot
kmerplotdat <- evalSeqs %>% filter(str_sub(Metric, start = -7) == "Entropy")
kmerplotdat$Metric <- gsub(" Entropy", "", kmerplotdat$Metric)
kmerplotdat$Metric <- factor(kmerplotdat$Metric, levels = c(
  "Sequence", "32mer", "24mer", "16mer", "8mer"))

p2b <- ggplot(kmerplotdat, aes(x=Metric, y=Value, color=Dataset, group=Dataset)) +
  geom_line() +
  #geom_point(size=1) +
  scale_color_manual(values=pal8) +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  labs(x="", y="Sequence Entropy")


## stitch two together plots with 2 col, 1row grid:
ggarrange(p2a, p2b, common.legend = TRUE, align = "h")
setwd("~/github/COIdatabases/rescript_paper_data/")   ## change as needed
ggsave("boldCOI_chordateOnly_SequenceEval_Figure.png", width=20, height=10, units="cm")

