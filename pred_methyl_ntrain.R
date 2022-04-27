# pred_model_ntrain.R
# -----------------------------------------------------------------------------
# Author:             Chang Chen
# Date last modified: April 1, 2022
#
# Run prediction pipeline

# Packages ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(janitor)
  library(supersigs)
})
source(here("./code/methyl_analysis_new_helpers.R"))

# find all the folders
# folders <- list.dirs(here('code/'),recursive=FALSE)

# Parameters ----
use_supersigs = FALSE # combine supersigs or not
use_clusters = TRUE
methylation = "enhancer" # "enhancer", "promoter", "both", or "all"
signature_dir = "UV*_SKCM"
type = str_split(signature_dir,"_")[[1]][1]
tissue = str_split(signature_dir,"_")[[1]][2]

# Read data -----
annoFolder <- here('data/Anno/')
file <- list.files(annoFolder,pattern=tissue)
filePath <- file.path(annoFolder,file)
# get object
methObj = readRDS(filePath)
# select regions
if(methylation=="enhancer"){
  methObj = methObj[!is.na(rowData(methObj)$value_enhancers_fantom), ]
}else if(methylation=="promoter"){
  methObj = methObj[!is.na(rowData(methObj)$value_promoters), ]
}else if(methylation=="both"){
  methObj = methObj[!is.na(rowData(methObj)$value_enhancers_fantom)|!is.na(rowData(methObj)$value_promoters), ]
}else if(methylation=="all"){
}else{
  stop("please input a methylation region enhancer, promoter, both, or all.")
}

if(use_clusters){
  p_clusters = readRDS(here(paste0("data/Clusters/p_clusters_",signature_dir,"_",methylation,".rds")))
  methObj <- clustered_se(methObj,p_clusters)
}

# Get signature_dt for
# (SMOKING - BLCA) = 31 (SMOKING - LUAD) = 32 (OBESITY - UCEC) = 55
# (OBESITY - KIRP) = 56 (OBESITY - ESCA) = 57 (OBESITY - COAD) = 58
# (ALCOHOL - LIHC) = 61 (POLE    - BRCA) = 42 (POLE    - UCEC) = 39
signature_caf = readRDS(here("./data/signature_caf.rds"))
signame <- colnames(signature_caf)[tissue == signature_caf[1,]& signature_caf[4,]==type]
if(tissue=="ESCA"&type=="SMOKING"){signame="SMOKING (ESCAD)"}
mutation_dt <- signature_caf[["Data",signame]][["DataSetFiltered"]]%>%
  dplyr::select('PATIENT','AGE','IndVar',contains(c('[',']','>')))%>%
  dplyr::rename(patient=PATIENT) %>%
  distinct()

mutation_dt$sample_id=1:nrow(mutation_dt)
colData(methObj) = colData(methObj) %>% 
  as_tibble() %>%
  left_join(., mutation_dt%>%dplyr::select(patient,IndVar)%>%clean_names(), by="patient") %>% DataFrame()
rownames(colData(methObj)) = colData(methObj)$patient

# remove patient with NA ind_var
methObj = methObj[,!is.na(methObj$ind_var)]

# 5 * 5-fold cross validations
library(doParallel)  
# no_cores <- detectCores() - 1  
cl <- makeCluster(5)
registerDoParallel(cl) 
set.seed(100)
out = list()
for(j in 1:5){ # 5 iterations
  model_dt = t(assay(methObj)) %>%
    as_tibble() %>%
    mutate(ind_var = factor(colData(methObj)$ind_var),
           patient=methObj$patient)
  folds <- createFolds(methObj$ind_var, k = 5) # 5-fold CV
  for(i in 1:5){
    message(paste0(j,"-",i))
    message("inner 5*3-fold-cross validation for selecting optimal n")
  # Step 1: inner 5*3-fold-cross-validation to find optimal n
    inner.out = list()
    for(s in 1:5){
      index <- setdiff(1:ncol(methObj), folds[[i]])
      inner_folds <- createFolds(methObj$ind_var[index], k = 3)
      for(k in 1:3){
        message(paste0("...testing CV ",s," fold ", k))
        train_index <- setdiff(index,inner_folds[[k]])
        results_limma = run_limma(methObj[,train_index]) 
        if(use_supersigs==FALSE){
          tryCatch(
            { 
              train = model_dt[train_index,]
              test = model_dt[-train_index,]
              results_pred <- suppressWarnings(pred_func(train = train, test = test, results_limma))
              if(!is.null(results_pred)){inner.out = c(inner.out,list(results_pred$auc))}
            },
            error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
          )
        }else{
          tryCatch(
            {
              # get mutational signatures from training set
              mutation_dt_train = mutation_dt%>%
                dplyr::filter(patient%in%model_dt$patient[train_index])%>%
                dplyr::select(!c(patient)) 
              supersig = get_signature(data = mutation_dt_train, factor = type)
              if(!is.null(supersig)){
                mut_sigs = unlist(supersig@Features)
                names(mut_sigs) = mut_sigs}
              else{
                mut_sigs = NULL
                }
              # define training set and test set
              mutation_dt1 = mutation_dt %>% dplyr::select(all_of(mut_sigs),patient)%>%
                mutate_if(is.numeric,scale) %>% as.data.table()
              model_dt = as.data.table(model_dt)
              train = mutation_dt1[model_dt[train_index,],on="patient"]%>%as.data.frame()
              test = mutation_dt1[model_dt[-train_index,],on="patient"]%>%as.data.frame()
              results_pred_mut <- suppressWarnings(pred_func_mut_methyl(train = train, test = test, results_limma,mut_sigs))
              if(!is.null(results_pred_mut)){inner.out = c(inner.out,list(results_pred_mut$auc))}
            },
            error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
          )
          }
      }
    }
    if(is.null(unlist(inner.out))){break}
    inner.result = tibble(n_keep_genes=rep(c(5,10,25,50,100,200,500),times=length(inner.out)),auc=unlist(inner.out))
    auc.median = inner.result %>% group_by(n_keep_genes) %>% summarise(auc_m=median(auc))
    n_pick = auc.median[which.max(auc.median$auc_m),"n_keep_genes"]
    
    # Step 2: train on the three folds and test on the remaining fold
    train_index <- setdiff(1:ncol(methObj), folds[[i]])
    results_limma = run_limma(methObj[,train_index]) 
    if(use_supersigs==FALSE){
      tryCatch(
        { 
          train = model_dt[train_index,]
          test = model_dt[-train_index,]
          results_pred <- suppressWarnings(pred_func_fixn(train = train, test = test, results_limma,n_pick))
          out = append(out,list(results_pred))
        },
      error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
      )
    }else{
      tryCatch(
        {
          # get mutational signatures from training set
          mutation_dt_train = mutation_dt%>%
            dplyr::filter(patient%in%model_dt$patient[train_index])%>%
            dplyr::select(!c(patient)) 
          supersig = get_signature(data = mutation_dt_train, factor = type)
          mut_sigs = unlist(supersig@Features)
          names(mut_sigs) = mut_sigs
          # define training set and test set
          mutation_dt1 = mutation_dt %>% dplyr::select(all_of(mut_sigs),patient)%>%
            mutate_if(is.numeric,scale) %>% as.data.table()
          model_dt = as.data.table(model_dt)
          train = mutation_dt1[model_dt[train_index,],on="patient"]%>%as.data.frame()
          test = mutation_dt1[model_dt[-train_index,],on="patient"]%>%as.data.frame()
          results_pred_mut <- suppressWarnings(pred_func_mut_methyl_fixn(train = train, test = test, results_limma,mut_sigs,n_pick))
          out = append(out,list(results_pred_mut))
        },
        error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
      )}
  }
}

stopCluster(cl)  

output = matrix(unlist(out),ncol=2,byrow=TRUE)
colnames(output) = c("n_keep_probes","AUC")
print(paste0(signature_dir,"_",methylation,"_SuperSigs_",use_supersigs))
# mean AUC
print(paste0("mean: ",mean(output[,"AUC"])))
# sd AUC
print(paste0("sd: ",sd(output[,"AUC"])) 
# mode of n_keep_probes
names(which.max(table(output[,1])))
write.csv(output,file=here(paste0("code/",signature_dir,"/cv_scaled_ntrain_",methylation,"_supersigs_",use_supersigs,".csv")))
