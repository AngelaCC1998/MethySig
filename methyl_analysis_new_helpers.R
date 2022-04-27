# methyl_analysis_new_helpers.R
# -----------------------------------------------------------------------------
# Author:             Chang Chen
# Date last modified: Apr 27, 2022
#
# Helper functions for methylation signature analysis


suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(ROCR)
  library(Biobase)
  library(limma)
  library(scater)
  library(lumi)
  library(data.table)
  library(DMRcate)
})

# density plot for Beta values
dplot <- function(methObj){
  dp <- assays(methObj)[[1]] %>% as.data.frame() %>% 
    pivot_longer(cols=starts_with("TCGA"),names_to="patient",values_to = "Beta") %>%
    left_join(colData(methObj)%>%as.data.frame()%>%dplyr::select(patient,ind_var),"patient") %>%
    ggplot(aes(x=Beta)) +
   #   geom_jitter(height=0,width=0.2) +
  #  theme_bw()
      geom_histogram(aes(y = ..density..)) +
      facet_wrap(vars(ind_var)) 
  return(dp)
}

# density plot for Beta mean values
dmplot <- function(methObj){
  exp_ind <- which(methObj$ind_var==TRUE)
  exposed_patients <- methObj$patient[exp_ind]
  unexposed_patients <- methObj$patient[-exp_ind]
  exposed_probe_mean <- rowMeans(assay(methObj)[,exposed_patients])
  unexposed_probe_mean <- rowMeans(assay(methObj)[,unexposed_patients])
  dmp <- data.frame(value_exposed=exposed_probe_mean,value_unexposed=unexposed_probe_mean) %>%
    pivot_longer(cols=c(value_exposed,value_unexposed),names_to="exposure",names_prefix="value_",values_to="value")%>%
    ggplot(aes(x=value)) +
    geom_histogram(aes(y = ..density..)) +
    facet_wrap(vars(exposure)) 
  return(dmp)
}

# Get correlation matrices separately for exposed and unexposed patients
get_cor_matrices = function(methObj){
  h_matrix = t(assay(methObj))
  
  exposed_patients = as_tibble(colData(methObj)) %>%
    filter(ind_var == TRUE) %>%
    pull(patient)
  
  unexposed_patients = as_tibble(colData(methObj)) %>%
    filter(ind_var == FALSE) %>%
    pull(patient)
  
  cor_exposed_matrix = list()
  cor_unexposed_matrix = list()
  for(chr_i in paste0("chr", 1:22)){
    chr_genes_ordered = as_tibble(rowRanges(methObj)) %>%
      filter(seqnames == chr_i) %>%
      arrange(start) %>%
      pull(probes)
    
    cor_exposed_matrix[[chr_i]] = cor(h_matrix[exposed_patients, chr_genes_ordered])
    cor_unexposed_matrix[[chr_i]] = cor(h_matrix[unexposed_patients, chr_genes_ordered])
  }
  
  out = list(cor_exposed_matrix, cor_unexposed_matrix)
  names(out) = c("Exposed", "Unexposed")
  return(out)
}

# Get hierarchical clusters from correlation matrices
get_hclusters = function(cors){
  cor_exposed_matrix = cors[["Exposed"]]
  cor_unexposed_matrix = cors[["Unexposed"]]
  
  row_clusters_exposed = list()
  row_clusters_unexposed = list()
  
  for(chr_i in paste0("chr", 1:22)){
    # Remove NAs
    nona_genes = rowSums(is.na(cor_exposed_matrix[[chr_i]])) != ncol(cor_exposed_matrix[[chr_i]]) - 1
    nona_genes = nona_genes & (rowSums(is.na(cor_unexposed_matrix[[chr_i]])) != ncol(cor_unexposed_matrix[[chr_i]]) - 1)
    
    # Exposed correlation matrix
    cor_exposed_nona_matrix = cor_exposed_matrix[[chr_i]][nona_genes, nona_genes]
    row_clusters_exposed[[chr_i]] = hclust(dist(cor_exposed_nona_matrix))
    
    # Unexposed correlation matrix
    cor_unexposed_nona_matrix = cor_unexposed_matrix[[chr_i]][nona_genes, nona_genes]
    row_clusters_unexposed[[chr_i]] = hclust(dist(cor_unexposed_nona_matrix))
  }
  
  out = list(row_clusters_exposed, row_clusters_unexposed)
  names(out) = c("Exposed", "Unexposed")
  return(out)
}

# Get clusters with high correlation from saved clusters dendrogram in heatmap
extract_clusters = function(clusters_input, cors_input){
  clusters_tb = NULL
  k = 1
  
  # Get clusters with high correlation and low sparsity
  for(chr_i in paste0("chr", 1:length(cors_input))){
    max_k = nrow(cors_input[[chr_i]])/10 # max number of clusters
    corr_summ = tibble(gene = rownames(cors_input[[chr_i]]),
                       mean_corr = apply(abs(cors_input[[chr_i]]), 1, function(x) mean(x, na.rm = TRUE)))
    
    for(i in 1:max_k){
      clusters = cutree(clusters_input[[chr_i]], k = i)
      for(j in 1:i){
        genes_in_j = names(clusters[clusters == j])
        if(!is.null(clusters_tb)){
          if(any(genes_in_j %in% clusters_tb$gene))
            next
        }
        cor_subset = cors_input[[chr_i]][genes_in_j, genes_in_j]
        avg_cor = mean(cor_subset)
        
        if(avg_cor > 0.5){
          cluster_tb = tibble(gene = genes_in_j,
                              cluster = k,
                              avg_cor = avg_cor)
          clusters_tb = bind_rows(clusters_tb, cluster_tb)
          k = k+1
        }
      }
    }
  }
  
  return(clusters_tb)
}


# Get partitioning of exposed and unexposed clusters
partition_clusters = function(exposed_clusters, unexposed_clusters){
  # If no clusters in either or both
  if(is.null(exposed_clusters) && is.null(unexposed_clusters)){
    return(NULL)
  } else if(is.null(exposed_clusters)){
    p_clusters = unexposed_clusters %>%
      mutate(partitioned_cluster = paste("cluster", cluster, "NA", sep = "-")) %>%
      dplyr::select(gene, partitioned_cluster)
    return(p_clusters)
  } else if(is.null(unexposed_clusters)){
    p_clusters = exposed_clusters %>%
      mutate(partitioned_cluster = paste("cluster", "NA", cluster, sep = "-")) %>%
      dplyr::select(gene, partitioned_cluster)
    return(p_clusters)
  } else {
    p_clusters = exposed_clusters %>%
      full_join(., unexposed_clusters, by = "gene") %>%
      mutate(partitioned_cluster = paste("cluster", cluster.x, cluster.y, sep = "-")) %>%
      dplyr::select(gene, partitioned_cluster)
    return(p_clusters)
  }
}


# get clustered mean values
clustered_se <- function(methObj,p_clusters){
  if(is.null(p_clusters)){
    return(methObj)
  }
  cD = colData(methObj)
  rD = rowData(methObj) %>% as.data.frame() %>%
    left_join(.,p_clusters,by="probes") %>%
    mutate(partitioned_cluster = ifelse(is.na(partitioned_cluster), probes, partitioned_cluster))
    
  as = data.frame(assay(methObj),check.names=FALSE) %>% 
    mutate(probes=rownames(methObj)) %>% 
    left_join(.,p_clusters,by="probes") %>%
    mutate(partitioned_cluster = ifelse(is.na(partitioned_cluster), probes, partitioned_cluster)) %>%  as.data.table()
  as2 = as[,-"probes"][,lapply(.SD,mean),by=partitioned_cluster] # calculate mean beta values for clusters for each patient
  aS = sapply(as2[,-"partitioned_cluster"],as.numeric)
  rownames(aS) = as2$partitioned_cluster
  # setkey(as,"partitioned_cluster")
  # setkey(as2,"partitioned_cluster")
  # as3 = as[,c("probes","partitioned_cluster")][as2]
  # aS = as3[,-c("probes","partitioned_cluster")] %>% as.data.frame()
  # rownames(aS) = as3$probes
  # se = SummarizedExperiment(assays=aS,colData=cD,rowData=rD)
  se = SummarizedExperiment(assays=aS,colData=cD)
  return(se)
}


# probe-level methylation analysis
run_limma <- function(methObj){
  # convert Beta values to M values
  assays(methObj)[[1]] <- beta2m(assays(methObj)[[1]])
  # create design matrix
  IndVar <- factor(tolower(colData(methObj)$ind_var))
  design <- model.matrix(~0+IndVar) 
  colnames(design) <- levels(IndVar)
  # fit the linear model for each gene given a series of arrays
  fit <- lmFit(assays(methObj)[[1]], design)
  # create a contrast matrix for specific comparisons
  contMatrix <- makeContrasts(false-true,levels=design)
  # fit the contrasts
  fit2 <- contrasts.fit(fit, contMatrix)
  fit2 <- eBayes(fit2)  #empirical bayes
  # look at the numbers of DM CpGs at FDR < 0.05
  # print(summary(decideTests(fit2)))
  if(ncol(rowData(methObj))==0){
    annsub = data.frame(probes=rownames(methObj))
    DMPs <- topTable(fit2,num=Inf,coef=1,sort.by="p",genelist=annsub)
  }else{  
    annsub <- rowData(methObj)[match(rownames(assays(methObj)[[1]]),rowData(methObj)$probes),c(1,3:8)]
    DMPs <- topTable(fit2,num=Inf,coef=1,sort.by="p",genelist=annsub)
  }
  return(DMPs=DMPs)
}

# region-level methylation analysis
run_dmrcate <- function(methObj){
  # convert Beta values to M values
  assays(methObj)[[1]] <- beta2m(assays(methObj)[[1]])
  # create design matrix
  IndVar <- factor(tolower(colData(methObj)$ind_var))
  design <- model.matrix(~0+IndVar) 
  colnames(design) <- levels(IndVar)
  # create a contrast matrix for specific comparisons
  contMatrix <- makeContrasts(false-true,levels=design)
  # fit the limma for individual CpG sites
  myAnnotation <- cpg.annotate(object = assay(methObj), datatype = "array", what = "M", 
                               analysis.type = "differential", design = design, 
                               contrasts = TRUE, cont.matrix = contMatrix, 
                               coef = "false - true", arraytype = "450K")
  # a kernel estimate against a null comparison to identify DMRs
  DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
  results.ranges <- extractRanges(DMRs)
  #ov <- findOverlaps(results.ranges,rowRanges(methObj_prom),type="any",select="all")
  #rowRanges(methObj)[ov@to,]$probes
  return(results.ranges=results.ranges)
}

# Predict on test set with different n_keep
pred_func = function(train, test, results_limma){
  auc_vec = c()
  n_keep_genes = sapply(c(5, 10, 25, 50, 100, 200, 500), function(x) min(x, nrow(results_limma)))
  for(n_keep in n_keep_genes){
    dm_probe = results_limma %>%
      # arrange(pvalue) %>%
      dplyr::slice(1:(min(n_keep, nrow(.))))
    
    # Subset to differential genes
    train_tmp = train[, c(dm_probe$probes, "ind_var")]
    test_tmp = test[, c(dm_probe$probes, "ind_var")]
    fitControl <- trainControl(method = "repeatedcv",
                               number = 2,
                               repeats = 5)
    plr <- train(ind_var ~ ., data = train_tmp,
                 method = "plr",
                 trControl = fitControl,
                 preProcess=c('center','scale','nzv'))
    # plr
    
    predictions = predict(plr, newdata = test_tmp, type = "prob")[, 2]
    pred_ROCR = prediction(predictions, test_tmp$ind_var)
    roc_ROCR <- performance(pred_ROCR, measure = "tpr", x.measure = "fpr")
    # plot(roc_ROCR, main = "ROC curve")
    # abline(a = 0, b = 1)
    auc_ROCR <- performance(pred_ROCR, measure = "auc")
    auc_vec = c(auc_vec, auc_ROCR@y.values[[1]])
  }
  
  auc_tb = tibble(n_keep_genes = n_keep_genes,
                  auc = auc_vec)
  
  return(auc_tb)
}

# Predict on test set with different n_keep, randomly select genes
pred_func_rand = function(train, test, results_limma){
  auc_vec = c()
  n_keep_genes = sapply(c(5, 10, 25, 50, 100, 200, 500), function(x) min(x, nrow(results_limma)))
  for(n_keep in n_keep_genes){
    dm_probe = results_limma %>%
      # arrange(pvalue) %>%
      dplyr::slice_sample(n=min(n_keep, nrow(.)))
    
    # Subset to differential genes
    train_tmp = train[, c(dm_probe$probes, "ind_var")]
    test_tmp = test[, c(dm_probe$probes, "ind_var")]
    fitControl <- trainControl(method = "repeatedcv",
                               number = 2,
                               repeats = 5)
    plr <- train(ind_var ~ ., data = train_tmp,
                 method = "plr",
                 trControl = fitControl,
                 preProcess=c('center','scale','nzv'))
    # plr
    
    predictions = predict(plr, newdata = test_tmp, type = "prob")[, 2]
    pred_ROCR = prediction(predictions, test_tmp$ind_var)
    roc_ROCR <- performance(pred_ROCR, measure = "tpr", x.measure = "fpr")
    # plot(roc_ROCR, main = "ROC curve")
    # abline(a = 0, b = 1)
    auc_ROCR <- performance(pred_ROCR, measure = "auc")
    auc_vec = c(auc_vec, auc_ROCR@y.values[[1]])
  }
  
  auc_tb = tibble(n_keep_genes = n_keep_genes,
                  auc = auc_vec)
  
  return(auc_tb)
}

# Predict on test set with different n_keep probes + mutational signatures
pred_func_mut_methyl = function(train, test, results_limma, mut_sigs){
  auc_vec = c()
  n_keep_genes = sapply(c(5, 10, 25, 50, 100, 200, 500), function(x) min(x, nrow(results_limma)))
  for(n_keep in n_keep_genes){
    dm_probe = results_limma %>%
      # arrange(pvalue) %>%
      dplyr::slice(1:(min(n_keep, nrow(.))))
    
    # Subset to differential genes
    train_tmp = train[, c(dm_probe$probes, mut_sigs, "ind_var")]
    test_tmp = test[, c(dm_probe$probes, mut_sigs, "ind_var")]
    fitControl <- trainControl(method = "repeatedcv",
                               number = 2,
                               repeats = 5)
    plr <- train(ind_var ~ ., data = train_tmp,
                 method = "plr",
                 trControl = fitControl,
                 preProcess=c('center','scale','nzv'))
    # plr
    
    predictions = predict(plr, newdata = test_tmp, type = "prob")[, 2]
    pred_ROCR = prediction(predictions, test_tmp$ind_var)
    roc_ROCR <- performance(pred_ROCR, measure = "tpr", x.measure = "fpr")
    # plot(roc_ROCR, main = "ROC curve")
    # abline(a = 0, b = 1)
    auc_ROCR <- performance(pred_ROCR, measure = "auc")
    auc_vec = c(auc_vec, auc_ROCR@y.values[[1]])
  }
  
  auc_tb = tibble(n_keep_genes = n_keep_genes,
                  auc = auc_vec)
  
  return(auc_tb)
}

# Predict on test set with fixed n_keep
pred_func_fixn = function(train, test, results_deseq, n_keep){
  dm_genes = results_deseq %>%
    #   arrange(pvalue) %>%
    dplyr::slice(1:(min(n_keep, nrow(.))))
  
  # Subset to differential genes
  train_tmp = train[, c(dm_genes$probes, "ind_var")]
  test_tmp = test[, c(dm_genes$probes, "ind_var")]
  fitControl <- trainControl(method = "repeatedcv",
                             number = 2,
                             repeats = 5)
  plr <- train(ind_var ~ ., data = train_tmp,
               method = "plr",
               trControl = fitControl,
               preProcess=c('center','scale','nzv'))
  
  predictions = predict(plr, newdata = test_tmp, type = "prob")[, 2]
  pred_ROCR = prediction(predictions, test_tmp$ind_var)
  roc_ROCR <- performance(pred_ROCR, measure = "tpr", x.measure = "fpr")
  # plot(roc_ROCR, main = "ROC curve")
  # abline(a = 0, b = 1)
  auc_ROCR <- performance(pred_ROCR, measure = "auc")
  
  auc_tb = tibble(n_keep_genes = n_keep,
                  auc = auc_ROCR@y.values[[1]])
  
  return(auc_tb)
}

# Predict on test set with fixed n_keep probes + mutational signatures
pred_func_mut_methyl_fixn = function(train, test, results_limma, mut_sigs,n_keep){
    dm_probe = results_limma %>%
      # arrange(pvalue) %>%
      dplyr::slice(1:(min(n_keep, nrow(.))))
    
    # Subset to differential genes
    train_tmp = train[, c(dm_probe$probes, mut_sigs, "ind_var")]
    test_tmp = test[, c(dm_probe$probes, mut_sigs, "ind_var")]
    fitControl <- trainControl(method = "repeatedcv",
                               number = 2,
                               repeats = 5)
    plr <- train(ind_var ~ ., data = train_tmp,
                 method = "plr",
                 trControl = fitControl,
                 preProcess=c('center','scale','nzv'))
    
    predictions = predict(plr, newdata = test_tmp, type = "prob")[, 2]
    pred_ROCR = prediction(predictions, test_tmp$ind_var)
    roc_ROCR <- performance(pred_ROCR, measure = "tpr", x.measure = "fpr")
    # plot(roc_ROCR, main = "ROC curve")
    # abline(a = 0, b = 1)
    auc_ROCR <- performance(pred_ROCR, measure = "auc")
  
  auc_tb = tibble(n_keep_genes = n_keep,
                  auc = auc_ROCR@y.values[[1]])
  
  return(auc_tb)
}
