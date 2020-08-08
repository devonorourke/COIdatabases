## script to create figure plotting frequency of homopolymers for gap-removed raw COI sequences

library(tidyverse)
library(scales)

###################################################
## import data plot
###################################################

df <- read_table(file="https://raw.githubusercontent.com/devonorourke/COIdatabases/master/data/homopolymer_motifs.table", col_names = FALSE)
colnames(df) <- c("NumberOfSequencesWithHomopolymer", "HomopolymerLength")

ggplot(df, aes(x=HomopolymerLength, 
               y=NumberOfSequencesWithHomopolymer, 
               label=NumberOfSequencesWithHomopolymer)) +
  geom_col() +
  geom_text(size=2.5, vjust=-1, hjust=0, angle=12.5) +
  scale_y_continuous(labels = comma) +
  theme_light() +
  theme(panel.border = element_blank()) +
  labs(x="Homopolymer length (bp)", 
       y="Sequences with at least 1 homopolymer")

setwd("/Users/devonorourke/github/COIdatabases/figures")  ## change path as needed for your machine
ggsave("homopolymer_plot.png", width = 20, height = 15, units = "cm")
