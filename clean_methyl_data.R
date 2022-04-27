# clean_methyl_data.R
# -----------------------------------------------------------------------------
# Adapted from Albert's code
# Date last modified: Dec 24, 2021
#
# Clean methylation data (apply after download_methyl_data.R)
library(pacman)
p_load(here, tidyverse, SummarizedExperiment, annotatr, tools, data.table, janitor,TxDb.Hsapiens.UCSC.hg38.knownGene)

fullDataFolder = here("data/AllProbes/")
normalsFolder = here("data/Normals/")
noXYFolder = here("data/NoXY/")
AnnoFolder = here("data/Anno/")


# Load SummarizedExperiment
filepath = list.files(noXYFolder, full.names = T, pattern="MESO")
load(filepath)

obj = file_path_sans_ext(basename(filepath))
methObj = get(obj)

# remove probes with at least one NA
methObj = subset(methObj,subset = (rowSums(is.na(assay(methObj))) == 0))

# remove FFPE samples
# methObj = methObj[,!methObj$is_ffpe]
# remove sample that is not from primary tumor but from normal tissues, 
# that's why we have duplicate records for one patient in most of the cases
# methObj = methObj[, !methObj$shortLetterCode%in%c("NB","NT","NBC","NEBV","NBM")]
# remove duplicated records for one patient
dup = methObj$patient[duplicated(methObj$patient)]
drop.records = sapply(dup,function(d) {
  SL = methObj$shortLetterCode[methObj$patient==d]
  if('TP' %in% SL ) {
    if(unique(SL)!='TP'){which(methObj$patient==d & methObj$shortLetterCode!='TP')}
    else{which(methObj$patient==d)[sample(1:length(SL),length(SL)-1)]}
    }
  else if('TM' %in% SL){which(methObj$patient==d & methObj$shortLetterCode!='TM')}
  else {which(methObj$patient==d & methObj$shortLetterCode==SL[sample(1:length(SL),1)])}
})
methObj = methObj[,-unique(unlist(drop.records))]
# check if there're remianing duplicates
sum(duplicated(methObj$patient))

# Annotate probes (rowData)
methy_dat = rowRanges(methObj)
annots = c("hg38_genes_promoters","hg38_enhancers_fantom","hg38_basicgenes")
annotations = build_annotations(genome="hg38",annotations = annots)
dm_annotated = annotate_regions(
    regions = methy_dat,
    annotations = annotations)

annotations_dt = tibble(probes = names(dm_annotated),
                        annot = sub("(hg38_genes_)|(hg38_)", "", dm_annotated$annot$type),
                        gene = dm_annotated$annot$symbol,
                        value = 1)
# For each probe, annotate whether it is in the promoter/enhancer/etc. region
# of any gene (columns "value_" = # of genes)
# and what the list of genes are (columns "gene_")
annotations_wide_dt = annotations_dt %>%
    pivot_wider(names_from = annot, values_from = c(gene, value),
                values_fn = list("gene" = list, "value" = length))

# Add back to rowData
rowData(methObj) = rowData(methObj) %>%
    as_tibble(rownames = "probes") %>%
    left_join(., annotations_wide_dt, by = "probes")

# Add clinical information (colData)
# signature_caf = readRDS("signature_caf.rds")
# signature_dt = signature_caf[["Data", "BMI (UCEC)"]]$DataSetFiltered %>%
#   dplyr::select(PATIENT, AGE, GENDER, IndVar) %>%
#   clean_names() %>%
#   as_tibble()
# colData(methObj) = colData(methObj) %>% 
#     as_tibble() %>%
#     # dplyr::select(!c(age, gender.y, ind_var)) %>%
#     left_join(., signature_dt, by="patient") %>% DataFrame()

# Save annotated object
saveRDS(methObj, file = file.path(AnnoFolder, paste0(obj, "_annotated.rds")))


## wrap up the function
clean_methyl <- function(methObj, signature_dt, annotations){
  # Annotate probes (rowData)
  # annots = c("hg38_genes_promoters","hg38_enhancers_fantom","hg38_basicgenes")
  # annotations = build_annotations(genome="hg38",annotations = annots)
  methy_dat = rowRanges(methObj)
  dm_annotated = annotate_regions(
    regions = methy_dat,
    annotations = annotations)
  annotations_dt = tibble(probes = names(dm_annotated),
                          annot = sub("(hg38_genes_)|(hg38_)", "", dm_annotated$annot$type),
                          gene = dm_annotated$annot$symbol,
                          value = 1)
  annotations_wide_dt = annotations_dt %>%
    pivot_wider(names_from = annot, values_from = c(gene, value),
                values_fn = list("gene" = list, "value" = length))
  
  # Add back to rowData
  rowData(methObj) = rowData(methObj) %>%
    as_tibble(rownames = "probes") %>%
    left_join(., annotations_wide_dt, by = "probes")
  
  # Add clinical information (colData)
  colData(methObj) = colData(methObj) %>% 
    as_tibble() %>%
    left_join(., signature_dt, by="patient") %>% DataFrame()
  
  # Save annotated object
  saveRDS(methObj, file = file.path(annoFolder, paste0(obj, "_annotated.rds")))
}

