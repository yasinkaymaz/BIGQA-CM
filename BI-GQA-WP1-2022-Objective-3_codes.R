library(Seurat)
library(tidyverse)
library(pheatmap)
source("BIGQA-WP1-functions.R")

#The link to the seurat object is here: https://drive.google.com/drive/folders/1yvK1RHT4fut96rFHQkHdrQGmqsxvThrw 
seuObj <- readRDS("Bachireddy.all.tcells.seurat.rds")

#Clusters to examine:
cls <- names(table(seuObj@meta.data$cluster_id))#These are all Biscuit clusters from the paper. Stores redundancy. 
# or optionally a small subset
cls <- c(1, 7, 14, 2)
#To create the clusters distance matrix, cdf, run this wrapper function:
cdm <- ClusterDistances(seu = seuObj, cls_list = cls)

#' To visually inspect the distances on a heatmap (any alternative to pheatmap can work)
pheatmap::pheatmap(cdm)

#' Once Meta-clusters are defined, the `meta.data` can be modified to add a `meta_cluster` column.
#' Then, the function `WeiTtest` can be run on "All" samples. 

cluster_ratioStats <- WeiTtest(Pid = seu$patient_response,
                               Groups = seu$Group,
                               Clusters = seu$meta_cluster,
                               nTry = 2000, N_r = 6, N_nr = 5,
                               GroupNames = c("R", "NR"))

#Example plot similar to Figure 2C - meta-cluster 1 (MC1: combination of 19th and 29th clusters. Note: There is a type in the manuscript. MC1 should have 29 instead of 28. Check Sup.Fig3C)
cluster_ratioStats@ratios %>% 
  filter(Cluster == 'MC1_19_29') %>% 
  ggplot(aes(Groups, koi, fill=Groups))+geom_boxplot()+
  ylab("Proportion of cells")

