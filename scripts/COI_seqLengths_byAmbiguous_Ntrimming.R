## script to create figure plotting frequency of sequence lengths of ...
## ... N-trimmed ('trimmed'), or not, ('untrimmed') COI sequences that had ...
## ... gaps removed

library(tidyverse)

trimmed <- read_table(col_names=FALSE, "https://raw.githubusercontent.com/devonorourke/COIdatabases/master/data/seqlength_Ntrimmed_hist.txt") %>% 
  mutate(Group = "trimmed")

untrimmed <- read_table(col_names = FALSE,
                        "https://raw.githubusercontent.com/devonorourke/COIdatabases/master/data/seqlength_untrimmed_hist.txt") %>% 
  mutate(Group = "untrimmed")

plotdat <- rbind(trimmed, untrimmed)
colnames(plotdat) <- c("Counts", "SeqLength", "Group")

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
