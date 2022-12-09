library(Seurat)
library(tidyverse)

####################### DEFINE THE FUNCTIONS  #######################

#' This function returns a list of two arrays. The first is the Pearson correlation coefficients of 
#' mean gene expression of genes in a given gene set (Gs). The second is the df of Factor loading values per cell.
#' GEM: gene expression matrix. Has to be genes (rows) by cells (cols).
#' Gs: Gene set. Gene symbols have to be present the GEM row names.
#' nF: Number of factors to be run
RunCFA_GScors <- function( GEM, Gs, nF=3){
  
  #1st part: Compute the common factors
  fct <- factanal(GEM, nF, rotation="varimax")
  factors <- as.data.frame(head(fct$loadings, nrow(GEM)))
  #2nd part: Compute the mean Gs expression per cell
  GS.Cell.mean <- apply(GEM[which(rownames(GEM) %in% Gs), ], 2, mean, na.rm=T)#this requires tidyverse. update.
  
  GS.F.cors <- c()
  for(i in 1:ncol(factors)){
    GS.F.cors[i] <- cor(GS.Cell.mean, factors[[i]])
  }
  return(list(FactCors=GS.F.cors, Factors=factors))
}

#' GEM: gene expression matrix. Has to be genes (rows) by cells (cols).
#' Obs.cor: Observed Pearson correlation coefficient between one Factor and mean gene expression of a GS
#' Gs: Gene set. Needed for only size
#' Factor: An array of factor loadings per cell corresponing to the GEM.
#' ntry: Number of permutation
permutationTest <- function(Obs.cor, Gs, GEM, Factor, ntry=500){
  set.seed(878712)
  nGns <- nrow(GEM)#number of genes.
  Gs.size <- length(Gs)
  #ncd: Null correlations distribution
  ncd <- c()
  #p_val <- 0
  for(i in 1:ntry){
    ncd[i] <- cor(apply(GEM[sample(nGns, Gs.size), ], 2, mean), Factor)
  }
  p <- sum(abs(ncd) >= abs(GS.F.cors[1]))/(ntry)
  if(p > 0.5){p_val <- 1-p}else{p_val <- p}
  results <- list(p_val=p_val, Ncd=ncd)
  return(results)
}


####################### CODE APPLICATION #######################

#Downloaded count matrix https://drive.google.com/drive/folders/1vMHAPbBSRouoyRymJ0cFdH1OBooYmK9f
#Cells were downsampled to 1000 to recapitulate the method in the Bachireddy et al.  
cnt_b1 <- read.csv("/Users/yasinkaymaz/Documents/BIGQAWP1/cmlnf_B1_dense.csv", header=T, row.names = 1)
cnt_b1 <- t(cnt_b1[sample(nrow(cnt_b1), 1000), ])
seub1 <- CreateSeuratObject(counts = cnt_b1, project = "B01", min.cells = 3, min.features = 200)
seub1@meta.data$Outcome <- c("Responder")

cnt_b3 <- read.csv("/Users/yasinkaymaz/Documents/BIGQAWP1/cmlnf_B3_dense.csv", header=T, row.names = 1)
cnt_b3 <- t(cnt_b3[sample(nrow(cnt_b3), 1000), ])
seub3 <- CreateSeuratObject(counts = cnt_b3, project = "B03", min.cells = 3, min.features = 200)
seub3@meta.data$Outcome <- c("Non-Responder")
head(seub3@meta.data)

#Standard data processing by Seurat
seu <- merge(seub1, y = seub3, add.cell.ids = c("B01", "B03"), project = "DLIstudy")
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

#Import the gene sets
gs1 <- read.delim("~/Downloads/GS.proinf.txt", header = F)
gs1$V1

#Compute the CFA and correlations with mean Gene set expressions
cfaCor.test <- RunCFA_GScors(GEM = seu[["RNA"]]@data[, sample(ncol(seu[["RNA"]]@data), 100)], Gs = gs1$V1)

#Modify Seurat object for further investigation
seu@meta.data <- cbind(seu@meta.data, cfaCor.test$Factors)
seu@meta.data %>% ggplot(aes(Factor1, Factor2, color=Outcome))+geom_point()

#Permutation test for observed correlation coefficients 
test.f1 <- permutationTest(Obs.cor = GS.F.cors[1], #For only first factor. Repeat for others.
                           GEM=seu[["RNA"]]@data, 
                           Factor=seu@meta.data$Factor1,
                           Gs = gs1$V1)
test.f1$p_val
hist(test.f1$Ncd)
abline(v=GS.F.cors[1], lwd=3, col="red")



