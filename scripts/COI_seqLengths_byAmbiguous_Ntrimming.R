## script to create figure plotting frequency of sequence lengths of ...
## ... N-trimmed ('trimmed'), or not, ('untrimmed') COI sequences that had ...
## ... gaps removed

library(tidyverse)

###################################################
## import data and reformat for plot
###################################################
trimmed <- read_table(col_names=FALSE, "https://raw.githubusercontent.com/devonorourke/COIdatabases/master/data/seqlength_Ntrimmed_hist.txt") %>% 
  mutate(Group = "trimmed")

untrimmed <- read_table(col_names = FALSE,
                        "https://raw.githubusercontent.com/devonorourke/COIdatabases/master/data/seqlength_untrimmed_hist.txt") %>% 
  mutate(Group = "untrimmed")

plotdat <- rbind(trimmed, untrimmed)
colnames(plotdat) <- c("Counts", "SeqLength", "Group")

###################################################
## plot and save
###################################################

ggplot(plotdat, aes(x=SeqLength, y=Counts, color=Group)) +
  geom_point() +
  geom_vline(xintercept = 658, linetype="dotted", color="gray25") +
  facet_wrap(~ Group, nrow = 2) +
  scale_color_manual(values=c("firebrick", "dodgerblue")) +
  scale_y_continuous(trans = "log10", 
                     breaks = c(10,100,1000,10000,100000,1000000),
                     labels = c("10", "100", "1,000", "10,000", "100,000", "1,000,000")) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x="\nCOI sequence length (bp)", y="Number of COI sequences\n")

setwd('/Users/devonorourke/github/COIdatabases/figures')  ## amend this as needed for your own machine
ggsave('COI_seqLengths_byAmbiguous_Ntrimming.png', height = 15, width=20, units = "cm")


###################################################
## calculate fraction remaining if trimming at defined length intervals
###################################################

seq_intevals <- seq(50, 8000, 50)

seqs_trimmed <- plotdat %>% filter(Group == "trimmed") %>% summarise(sum(Counts)) %>% pull()
seqs_untrimmed <- plotdat %>% filter(Group == "untrimmed") %>% summarise(sum(Counts)) %>% pull()

SeqsRetained_bySeqLengthBin_trimmed <- plotdat %>% 
  filter(Group == "trimmed") %>%
  mutate(Bins = cut_interval(x=SeqLength, length = 50)) %>% 
  mutate(Bins = gsub("\\(", "", Bins)) %>% 
  mutate(Bins = gsub("\\[", "", Bins)) %>% 
  mutate(Bins = gsub("\\]", "", Bins)) %>% 
  mutate(Bins = gsub(",.*$", "", Bins)) %>% 
  mutate(Bins = as.integer(Bins)) %>% 
  group_by(Bins) %>% 
  summarise(SeqCountsByBin = sum(Counts)) %>% 
  arrange(Bins) %>% 
  mutate(cumSeqCountsByBin = cumsum(SeqCountsByBin)) %>% 
  mutate(trim_fracBinCounts = 1 - (cumSeqCountsByBin / seqs_trimmed)) %>% 
  mutate(trim_fracBinCounts = round(trim_fracBinCounts, 3)) %>% 
  select(Bins, trim_fracBinCounts)

SeqsRetained_bySeqLengthBin_untrimmed <- plotdat %>% 
  filter(Group == "untrimmed") %>%
  mutate(Bins = cut_interval(x=SeqLength, length = 50)) %>% 
  mutate(Bins = gsub("\\(", "", Bins)) %>% 
  mutate(Bins = gsub("\\[", "", Bins)) %>% 
  mutate(Bins = gsub("\\]", "", Bins)) %>% 
  mutate(Bins = gsub(",.*$", "", Bins)) %>% 
  mutate(Bins = as.integer(Bins)) %>% 
  group_by(Bins) %>% 
  summarise(SeqCountsByBin = sum(Counts)) %>% 
  arrange(Bins) %>% 
  mutate(cumSeqCountsByBin = cumsum(SeqCountsByBin)) %>% 
  mutate(untrim_fracBinCounts = 1 - (cumSeqCountsByBin / seqs_untrimmed)) %>% 
  mutate(untrim_fracBinCounts = round(untrim_fracBinCounts, 3)) %>% 
  select(Bins, untrim_fracBinCounts)

## this table was used to generate the information presented in the RESCRIPt document
SeqsRetained_bySeqLengthBin_all <- merge(SeqsRetained_bySeqLengthBin_trimmed, 
                                         SeqsRetained_bySeqLengthBin_untrimmed)

## how many seqs remain in either trimmed or untrimmed if we filter at 200 bp?
plotdat %>% 
  filter(Group == "untrimmed") %>%
  filter(SeqLength >= 200) %>% summarise(TotalSeqs = sum(Counts))

plotdat %>% 
  filter(Group == "trimmed") %>%
  filter(SeqLength >= 200) %>% summarise(TotalSeqs = sum(Counts))
