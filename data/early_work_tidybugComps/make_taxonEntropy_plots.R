library(tidyverse)
library(scales)
library(scico)

setwd("~/Desktop/coi_qiime_tmp/")

readfunction <- function(filepathe){
  tmp_df <- read.delim(filepathe, skip = 2, header=FALSE, stringsAsFactors = FALSE)
  tmp_df <- tmp_df[-1]
  colnames(tmp_df) <- c("Rank", "nUnique", "TaxonEntropy",
                        "nFeaturesAtDepth", "pFeaturesAtDepth", 
                        "nFeaturesClassified", "pFeaturesClassified",
                        "nFeaturesUnclassified", "pFeaturesUnclassified", "Dataset")
  tmp_df
}

c97_df <- readfunction("tidybug_cluster_97.tsv")
c98_df <- readfunction("tidybug_cluster_98.tsv")
c99_df <- readfunction("tidybug_cluster_99.tsv")
c100_df <- readfunction("tidybug_cluster_100.tsv")

all_clustData <- rbind(c97_df, c98_df, c99_df, c100_df)
rm(c97_df, c98_df, c99_df, c100_df)

# write.table(all_clustData, 
#             file="tidybug_evaluate_taxonomy_data.tsv",
#             quote = FALSE, row.names = FALSE, col.names = TRUE)


################
## plot:

ggplot(all_clustData,
       aes(x=Rank, y=TaxonEntropy, color=Dataset)) +
  geom_line() +
  theme_bw() +
  scale_color_viridis_d()
  #geom_col() +
  #facet_wrap(~Dataset)

plotdat <- all_clustData %>% 
  pivot_longer(c(-Rank, -Dataset), names_to="Metric", values_to="Values") %>% 
  filter(Metric != "TaxonEntropy")

plotdat$Dataset <- factor(plotdat$Dataset, levels = c(
  "tidybug_c97", "tidybug_c98", "tidybug_c99", "tidybug_c100"
))

plotdat$Metric <- factor(plotdat$Metric, levels=c(
  "nFeaturesAtDepth", "pFeaturesAtDepth",
  "nFeaturesClassified", "pFeaturesClassified",
  "nFeaturesUnclassified", "pFeaturesUnclassified",
  "nUnique"
))

ggplot(plotdat, 
       aes(x=Rank, y=Values, fill=Dataset)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~Metric, scales = "free", ncol=4) +
  theme_bw() +
  guides(fill=guide_legend(nrow=4,byrow=TRUE)) +
  theme(legend.justification = c(1,0),
        legend.position = c(1,0),
        legend.box.margin = margin(c(0,0,0,0)),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(colour = NA, fill = NA),
        strip.text = element_text(size=10)) +
  scale_fill_viridis_d(direction = -1) +
  labs(x="", y="", fill="")
## save as 'plot_clustering_FeaturesComps'
