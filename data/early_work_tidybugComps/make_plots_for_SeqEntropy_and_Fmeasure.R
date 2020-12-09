library(tidyverse)
library(scales)
library(scico)
library(formattable)

dat_df <- read.csv(file="~/Desktop/coi_qiime_tmp/topOrder_data.csv")
count_df <- read.csv(file="~/Desktop/coi_qiime_tmp/top16order_seqCounts.csv")
dat_df <- merge(dat_df, count_df)
dat_df$logVals <- log10(dat_df$nSeqs)
dat_df$PlotLabel <- paste0(dat_df$Dataset, "n = ", dat_df$logVals)

ggplot(dat_df %>% filter(Metric == 'F.Measure'), 
       aes(x=Level, y=Value, group=Dataset, color=logVals)) +
  facet_wrap(~Dataset) +
  geom_point() +
  geom_line() +
  scale_colour_scico(palette = 'batlow', begin = 0.2,
                     breaks=c(4,4.60206,5.20412,5.80618), 
                     labels=c('10,000','40,000','160,000','640,000')) +
  scale_x_continuous(labels=c("Phylum", "Class", "Order", "Family", "Genus", "Species")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
        panel.grid.minor = element_blank()) +
  labs(x="\nRank level", y="F-measure\n", color="num Seqs")


###################

entropy_df <- read.delim("~/Desktop/coi_qiime_tmp/all_SeqViz.tsv",
                         header = FALSE)
colnames(entropy_df) <- c('Metric', 'Value', 'Taxa')

EntropyVals <- c('8mer Entropy', '16mer Entropy', '24mer Entropy', '32mer Entropy', 'Sequence Entropy')
# LengthVals <- c('Length 1%', 'Length 25%', 'Length 75%', 'Length 99%',
#                 'Length max', 'Length median', 'Length min')

entropy_plotdat <- entropy_df %>% filter(Metric %in% EntropyVals)
entropy_plotdat$Metric <- factor(entropy_plotdat$Metric,
                                 levels = EntropyVals)

## can plot so we facet each Order...
ggplot(entropy_plotdat, 
       aes(x=Metric, y=Value, fill=Metric)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~Taxa) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
        legend.position = "top") +
  labs(x="", y="Shannon's entropy\n", fill="") +
  scale_fill_scico_d(palette = "hawaii")

## or plot so we facet by kmers:
ggplot(entropy_plotdat, 
       aes(x=Taxa, y=Value)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~Metric, nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
        legend.position = "top") +
  labs(x="", y="Shannon's entropy\n", fill="")

## which Orders have the biggest gap between 8mer and sequence entropy?
diff_df <- entropy_plotdat %>% 
  filter(Metric == '8mer Entropy' | Metric == 'Sequence Entropy') %>% 
  pivot_wider(values_from = "Value", names_from = "Metric") %>%
  mutate(EntropyDiff = `Sequence Entropy` - `8mer Entropy`) %>% 
  arrange(-EntropyDiff)
  
## color the diff's as heatmap style
formattable(diff_df,list(
  EntropyDiff = color_tile("lightblue", "orange")
))