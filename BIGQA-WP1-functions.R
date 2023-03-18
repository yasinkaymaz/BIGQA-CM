
####################### DEFINE THE FUNCTIONS  #######################

#' This function returns a list of two arrays. The first is the Pearson correlation coefficients of 
#' mean gene expression of genes in a given gene set (Gs). The second is the df of Factor loading values per cell.
#' @param GEM gene expression matrix. Has to be genes (rows) by cells (cols).
#' @param Gs Gene set. Gene symbols have to be present the GEM row names.
#' @param nF Number of factors to be run
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

#' @param GEM gene expression matrix. Has to be genes (rows) by cells (cols).
#' @param Obs.cor Observed Pearson correlation coefficient between one Factor and mean gene expression of a GS
#' @param Gs Gene set. Needed for only size
#' @param Factor An array of factor loadings per cell corresponing to the GEM.
#' @param ntry Number of permutation
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
  p <- sum(ncd >= Obs.cor[1])/(ntry)
  if(p > 0.5){p_val <- 1-p}else{p_val <- p}
  results <- list(p_val=p_val, Ncd=ncd)
  return(results)
}


#' Weighted t-test for cell ratios of patients from two groups.
#' This function utilizes R library 'weights'
#' @param Pid patient ids vector for each data point, i.e. cell.
#' @param Groups groups vector for each data point, i.e. cell.
#' @param Clusters cluster ids vector for each data point, i.e. cell.
#' @param nTry number of re-sampling for bootstrapping
#' @param N_r number of patient samples in one group
#' @param N_nr number of patient samples in the second group
#' @param GroupNames
WeiTtest <- function(Pid, Groups, Clusters, nTry=3000, N_r, N_nr, GroupNames){
  
  seu.meta.data <- data.frame(Pid= Pid, Groups = Groups, Clusters = Clusters)
  
  seu.meta.data %>%
    group_by(Clusters) %>%
    mutate(Uk=n()) %>% #Uk: total number of T cells in the cluster k
    group_by(Groups) %>% 
    mutate(Np=n()) %>% #Np: total number of T cells in group
    group_by(Pid) %>% 
    mutate(Ni=n()) %>% #Ni: total number of T cells of patient i
    group_by(Clusters, Pid) %>% 
    mutate(Nk=n()) %>% #Nk: total number of T cells of patient i in cluster k
    mutate(koi=Nk/Ni, iop=Ni/Np) %>% #koi: proportion of cells, iop: T cell contribution to the group pool by patient i
    mutate(Wi=ifelse(Groups==GroupNames[1], N_r*iop, N_nr*iop)) %>%
    unique() %>% as.data.frame() -> meta.ttest
  
  cls.exclude <- droplevels(as.data.frame(table(meta.ttest$Groups, meta.ttest$Clusters))$Var2[as.data.frame(table(meta.ttest$Groups, meta.ttest$Clusters))$Freq == 0])
  cls_list <- unique(Clusters)[!unique(Clusters) %in% cls.exclude]
  cls_stats <- c()
  examined.cls.list <- c()
  st.list <- list()
  rt.df <- data.frame()
  sm.plist <- list()
  for(k in cls_list){
    
    A <- meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[1]), "koi"]
    B <- meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[2]), "koi"]
    wA <- meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[1]), "Wi"]
    wB <- meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[2]), "Wi"]
    
    if(length(A) == 1 & length(B) == 1 | length(A) == 1){print(k); cls.exclude <- c(cls.exclude, k); next;}
    
    wtd.t.test(
      A,
      B,
      weight = wA,
      weighty = wB,
      samedata = F,
      alternative = "greater") -> tt
    cls_p_val <- tt$coefficients['p.value']
    
    U_k <- unique(meta.ttest[meta.ttest$Clusters == k,]$Uk)#The size of the cluster k
    pdist <- c()
    for(b in 1:nTry){
      rand.meta <- seu.meta.data
      rand.meta$Clusters <- as.vector(rand.meta$Clusters)
      #Randomly select U_k number of cells from the pool of all samples and compute the p-value using the same weighted t.test.
      cls = 'rand'
      rand.meta[sample(1:nrow(rand.meta), size = U_k, replace = T), 'Clusters'] <- cls
      
      rand.meta %>%
        group_by(Clusters) %>%
        mutate(Uk=n()) %>% #Uk: total number of T cells in the cluster k
        group_by(Groups) %>% 
        mutate(Np=n()) %>% #Np: total number of T cells in group
        group_by(Pid) %>% 
        mutate(Ni=n()) %>% #Ni: total number of T cells of patient i
        group_by(Clusters, Pid) %>% 
        mutate(Nk=n()) %>% #Nk: total number of T cells of patient i in cluster k
        mutate(koi=Nk/Ni, iop=Ni/Np) %>% #koi: proportion of cells, iop: T cell contribution to the group pool by patient i
        mutate(Wi=ifelse(Groups=="Responder", N_r*iop, N_nr*iop)) %>%
        unique() %>% 
        filter(Clusters == cls) %>%
        as.data.frame() -> rand.meta.ttest
      
      a <- rand.meta.ttest[which(rand.meta.ttest$Clusters == cls & rand.meta.ttest$Groups == GroupNames[1]), "koi"]
      b <- rand.meta.ttest[which(rand.meta.ttest$Clusters == cls & rand.meta.ttest$Groups == GroupNames[2]), "koi"]
      wa <- rand.meta.ttest[which(rand.meta.ttest$Clusters == cls & rand.meta.ttest$Groups == GroupNames[1]), "Wi"]
      wb <- rand.meta.ttest[which(rand.meta.ttest$Clusters == cls & rand.meta.ttest$Groups == GroupNames[2]), "Wi"]
      
      if(length(a) == 1 & length(b) == 1){next}
      
      #Compute the p-value on randomly selected cells.
      try(
        wtd.t.test(
          a,
          b,
          weight = wa,
          weighty = wb,
          samedata = F,
          alternative = "greater") -> tt
      )
      pdist <- c(pdist, tt$coefficients['p.value'])
      
    }
    
    cls_Q_val <- table(pdist < cls_p_val)['TRUE']/nTry
    names(cls_Q_val) <- "class_Q_val"
    
    print(paste("cluster:", k, "pval:", cls_p_val, "Qval:", cls_Q_val))
    
    sm.plist[[k]] <- pdist
    st.list[k] <- cls_Q_val
    
    examined.cls.list <- c(examined.cls.list, k)
    
    df <- data.frame(cbind(Cluster=k,
                           rbind(meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[1]), ],
                                 meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[2]), ])
    ))
    rt.df <- rbind(rt.df, df)
  }
  
  rt.df$Groups <- factor(rt.df$Groups, c(GroupNames[1], GroupNames[2]))
  
  object <- new(Class = "ClusCompStats",
                stats = st.list,
                ratios= rt.df,
                plist = sm.plist)
  
  return(object)
}

#' A simple Bhattacharyya Distance metric function based on discrete probability distibutions.
#' @param xA a vector of probabilities for each bin
#' @param xB another vector of probabilities for each bin
bhattDist <- function(xA, xB){
  bc <- 0
  for(n in 1:length(xA)){
    bc <- bc+sqrt(xA[n]*xB[n])
  }
  if(bc < 1.0){
    bd <- sqrt(1-bc)
  }else{
    bd <- 0
  }
  return(bd)
}


#' This function is a wrapper to create cluster-to-clusters distance matrix using Bhattacharyya Distance metric.
#' @param seu This is a Seurat object that stores Biscuit normalized expression matrix in its `counts` slot. Important: its `meta.data` slot has to store cluster ids under `cluster_id` column.
#' @param cls_list The list of clusters to compare to each other.
ClusterDistances <- function(seu, cls_list){
  
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
  
  exp <- seu@assays$RNA@counts #Biscuit normalized expression matrix
  cls_out <- c()
  df <- data.frame()
  for(i in cls_list){
    for(j in cls_list[!cls_list %in% cls_out]){
      
      cls_pairs <- c(i,j)
      cls_a <- seu@meta.data %>% filter(cluster_id == cls_pairs[1]) %>% rownames() #cluster a cell ids
      cls_b <- seu@meta.data %>% filter(cluster_id == cls_pairs[2]) %>% rownames() #cluster b cell ids
      
      if(!c(is_empty(cls_a)|is_empty(cls_b)) ){
        
        cls_a_ls <- apply(exp, 1, function(x) bhattDist(xA = hist(x[cls_a], breaks=hist(x[c(cls_a, cls_b)], plot=F)$breaks, plot=F)$counts/length(cls_a),
                                                        xB = hist(x[cls_b], breaks=hist(x[c(cls_a, cls_b)], plot=F)$breaks, plot=F)$counts/length(cls_b)
                                                        )
                          )
        
        bd <- log(sum(cls_a_ls))
        
      }else{
        bd <- NA
        print("not enough cells in the cluster")
      }
      
      df <- rbind(df, c(clsA=i, clsB=j, bd=bd))
      
    }
    cls_out <- c(cls_out, i)
  }
  
  colnames(df) <- c("A", "B", "BD")
  df[is.na(df)] <- 0
  cdf <- reshape::cast(df, A~B,value = "BD") 
  rownames(cdf) <- cdf[, 1]
  cdf <- cdf[, -c(1)]
  
  return(cdf)
}
