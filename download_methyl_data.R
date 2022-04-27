# download_methyl_data.R
# -----------------------------------------------------------------------------
# Author:             Chang Chen
# Date last modified: Dec 9, 2021
#
# Download methylation data (modified from Luda's code: downloadMethData_for_sharing.r)
# library(usethis) 
# usethis::edit_r_environ()
# Packages ---
library(TCGAbiolinks)
library(SummarizedExperiment)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(tidyverse)
library(data.table)
library(here)
source(here("code/downloadDataFunctions.r"))

fullDataFolder = here("data/AllProbes/")
normalsFolder = here("data/Normals/")
noXYFolder = here("data/NoXY/")

if(!file.exists(fullDataFolder)) dir.create(fullDataFolder)
if(!file.exists(normalsFolder)) dir.create(normalsFolder)
if(!file.exists(noXYFolder)) dir.create(noXYFolder)

#===================================
# list of cancer types to download
# cancerTypes = c('LAML','PAAD','ESCA','KIRP','LUAD','BLCA','BRCA','CESC','COAD','HNSC','KIRC','LGG','LIHC','PRAD','SKCM','LUSC','STAD','THCA','UCEC')
cancerTypes = c('LGG',"MESO")
# download and save hg38 methylation from harmonized data portal into fullDataFolder
sapply(cancerTypes, download450k, fullDataFolder)

#================================
# find probes to remove (X,Y chromosome probes, SNPs)
probesToExcl = getProbesToExclude()
# get paths to meth objects
methObjs = list.files(fullDataFolder, full.names = T, pattern="MESO")
# remove XY probes and save normals separately
sapply(methObjs, removeProbesAndSaveNorm, normFolder = normalsFolder, outputFolder = noXYFolder, probesToRemove = probesToExcl)

#====================================
#  create objects with no NA
# methObjs = list.files(noXYFolder, full.names = T)
# # remove XY probes and save normals separately
# sapply(methObjs[2:19],removeNAprobes,outputFolder = 'NoXY_noNA/')