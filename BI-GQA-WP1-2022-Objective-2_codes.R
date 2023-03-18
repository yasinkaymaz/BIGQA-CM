library(tidyverse)
library(Seurat)
library(weights)

source("BIGQA-WP1-functions.R")

#The link to the seurat object is here: https://drive.google.com/drive/folders/1yvK1RHT4fut96rFHQkHdrQGmqsxvThrw 
seuObj <- readRDS("Bachireddy.all.tcells.seurat.rds")

#Figure 2B is for only Pre-DLI samples. Therefore, need to subset.
Pre_samples <- c("B1", "B3", "B5", "B7", "B9", "B17", "B19", "B21", "B23", "B25", "B27")
seu.sub <- seuObj@meta.data %>% filter(sample_id %in% Pre_samples)

cluster_ratioStats <- WeiTtest(Pid = seu.sub$patient_response,
                               Groups = seu.sub$Group,
                               Clusters = seu.sub$cluster_id,
                               nTry = 2000, N_r = 6, N_nr = 5,
                               GroupNames = c("R", "NR"))

#Example plot similar to Figure 2B - cluster 4
cluster_ratioStats@ratios %>% 
  filter(Cluster == 4) %>% 
  ggplot(aes(Groups, koi, fill=Groups))+geom_boxplot()+
  ylab("Proportion of cells Pre-DLI")