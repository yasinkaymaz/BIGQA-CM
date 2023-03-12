
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
#' @param nTry number of re-sampling for boostrapping
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
  print(table(meta.ttest$Groups, meta.ttest$Clusters))
  print(cls.exclude)
  print(cls_list)
  cls_stats <- c()
  for(k in cls_list){
    
    a <- meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[1]), "koi"]
    b <- meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[2]), "koi"]
    wa <- meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[1]), "Wi"]
    wb <- meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[2]), "Wi"]
    
    wtd.t.test(
      a,
      b,
      weight = wa,
      weighty = wb,
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
      rand.meta[sample(1:nrow(rand.meta), size = U_k), 'Clusters'] <- cls
      
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
      #Compute the p-value on randomly selected cells.
      wtd.t.test(
        rand.meta.ttest[which(rand.meta.ttest$Clusters == cls & rand.meta.ttest$Groups == GroupNames[1]), "koi"],
        rand.meta.ttest[which(rand.meta.ttest$Clusters == cls & rand.meta.ttest$Groups == GroupNames[2]), "koi"],
        weight = rand.meta.ttest[which(rand.meta.ttest$Clusters == cls & rand.meta.ttest$Groups == GroupNames[1]), "Wi"],
        weighty = rand.meta.ttest[which(rand.meta.ttest$Clusters == cls & rand.meta.ttest$Groups == GroupNames[2]), "Wi"],
        samedata = F,
        alternative = "greater") -> tt
      
      pdist <- c(pdist, tt$coefficients['p.value'])
      
    }
    #hist(pdist)
    cls_Q_val <- table(pdist < cls_p_val)['TRUE']/nTry
    names(cls_Q_val) <- "class_Q_val"
    
    cls_stats <- c(cls_stats, cls_Q_val)
  }
  names(cls_stats) <- cls_list
  return(cls_stats)
}
