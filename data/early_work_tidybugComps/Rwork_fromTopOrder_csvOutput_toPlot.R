-modify format of file for plot 

##### start here....

library(stringr)
library(ggplot2)
library(tidyr)

setwd("/scratch/dro49/qiimetmp/tidybug_crossval")
orderlist <- read.csv("topOrder.list", header=FALSE, stringsAsFactors=FALSE)
taxalist <- str_replace(orderlist$V1, "o__", "")
filenames <- paste(getwd(),"/", "tidybug_tmpdir_",taxalist,"_evalClassifications/data.tsv",sep="")

dat_list <- lapply(filenames,function(i){
  tmp <- read.delim(i, header=TRUE, skip=2)
  tmp <- tmp[,-1]
  colnames(tmp) <- c("Level", "Precision", "Recall", "F.Measure", "Dataset")
  taxaname_1 <- str_replace(i, "/scratch/dro49/qiimetmp/tidybug_crossval/tidybug_tmpdir_", "")
  taxaname_2 <- str_replace(taxaname_1, "_evalClassifications/data.tsv", "")
  tmp$Dataset <- taxaname_2
  tmp
})

dat_df <- do.call("rbind", dat_list) %>%
pivot_longer(c(-Dataset, -Level), names_to="Metric", values_to="Value")

rm(dat_list, taxalist, orderlist, filenames)

write.csv(dat_df, file="topOrder_data.csv", quote=FALSE, row.names = FALSE)


########## to plot:
library(tidyverse)
library(scales)
library(scico)

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
