## plot used to generate bar chart of sequence lengths of COI fragments ...
## ... after primer-based coordinate trimming of giant_alignment file

library(ggplot2)
library(dplyr)
library(scales)

################################################################################
## load data and plot
################################################################################

df <- read_delim(file="https://raw.githubusercontent.com/devonorourke/COIdatabases/master/data/anml_seq_lengths.freq.table", 
                 col_names=FALSE, delim=" ")
colnames(df) <- c("NumberOfSeqs", "SeqLength")
df$SeqLength <- as.integer(df$SeqLength)
df$NumberOfSeqs <- as.integer(df$NumberOfSeqs)
sumSeqs <- sum(df$NumberOfSeqs)
df <- df %>%
  mutate(SeqLengtf = as.numeric(SeqLength)) %>%
  arrange(SeqLength) %>%
  mutate(CumSum = cumsum(NumberOfSeqs)) %>%
  mutate(fracSeqs = CumSum/sumSeqs)

ggplot(df, aes(x=SeqLength, y=NumberOfSeqs)) +
  geom_bar(stat="identity") +
  theme_bw() +
  scale_y_continuous(labels=comma) +
  labs(x="\nSequence length (bp)", y="Number of sequences\n")

setwd('/Users/devonorourke/github/COIdatabases/figures')
ggsave("anml_primer_COI_seqLength.png", height=10, width=15, units="cm")
