library(Seurat)
library(tidyverse)
source("~/Dropbox/codes/BIGQA-CM/BIGQA-WP1-functions.R")

####################### CODE APPLICATION #######################

#Downloaded count matrix https://drive.google.com/drive/folders/1vMHAPbBSRouoyRymJ0cFdH1OBooYmK9f
#Cells were downsampled to 1000 to recapitulate the method in the Bachireddy et al.  
cnt_b1 <- read.csv("/Users/yasinkaymaz/Documents/BIGQAWP1/cmlnf_B1_dense.csv", header=T, row.names = 1)
cnt_b1 <- t(cnt_b1[sample(nrow(cnt_b1), 1000), ])
seub1 <- CreateSeuratObject(counts = cnt_b1, project = "B01", min.cells = 3, min.features = 200)
seub1@meta.data$Outcome <- c("Responder")
rm(cnt_b1)
gc()

cnt_b3 <- read.csv("/Users/yasinkaymaz/Documents/BIGQAWP1/cmlnf_B3_dense.csv", header=T, row.names = 1)
cnt_b3 <- t(cnt_b3[sample(nrow(cnt_b3), 1000), ])
seub3 <- CreateSeuratObject(counts = cnt_b3, project = "B03", min.cells = 3, min.features = 200)
seub3@meta.data$Outcome <- c("Non-Responder")
head(seub3@meta.data)
rm(cnt_b3)
gc()

#Standard data processing by Seurat
seu <- merge(seub1, y = seub3, add.cell.ids = c("B01", "B03"), project = "DLIstudy")
rm(seub1, seub3);gc()
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:10)
DimPlot(seu, reduction = "umap", group.by = "Outcome")
DimPlot(seu, reduction = "pca", group.by = "Outcome")
rm(all.genes);gc()
#Import the gene sets
gs1 <- read.delim("GS.proinf.txt", header = F)
gs1$V1

#Compute the CFA and correlations with mean Gene set expressions
cfaCor.test <- RunCFA_GScors(GEM = seu[["RNA"]]@data[, sample(ncol(seu[["RNA"]]@data), )], Gs = gs1$V1)
saveRDS(cfaCor.test, file = "~/Documents/BIGQAWP1/cfaCor.test.rds")
#Modify Seurat object for further investigation
seu@meta.data <- cbind(seu@meta.data, cfaCor.test$Factors)
seu@meta.data %>% ggplot(aes(Factor1, Factor2, color=Outcome))+geom_point()

#Permutation test for observed correlation coefficients 
test.f1 <- permutationTest(Obs.cor = cfaCor.test[1]$FactCors[1], #For only first factor. Repeat for others.
                           GEM=seu[["RNA"]]@data, 
                           Factor=seu@meta.data$Factor1,
                           Gs = gs1$V1)
test.f1$p_val
hist(test.f1$Ncd)
abline(v=cfaCor.test[1]$FactCors[1], lwd=3, col="red")



