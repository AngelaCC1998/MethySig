# function to download TCGA data using TCGAbiolinks

geti = function(x,i)return(x[i])
# downloads 450k methylation from harmonized data portal (hg38)
# and saves in the specified folder
# if barcodes are specified, downloads only specified samples
download450k = function(cancer, folder = '', barcodes = NULL, removeFiles = T)
{
	library(TCGAbiolinks)
	library(SummarizedExperiment)
	# make a query
	# save object
	tryCatch({
		if(is.null(barcodes)) 
		{
			query_met.hg38 <- GDCquery(project= paste0("TCGA-",cancer), 
			                           data.type="Methylation Beta Value",
								   data.category = "DNA Methylation", 
								   platform = "Illumina Human Methylation 450")
		}else{
			query_met.hg38 <- GDCquery(project= paste0("TCGA-",cancer), 
								   data.category = "DNA Methylation", 
								   data.type="Methylation Beta Value",
								   platform = "Illumina Human Methylation 450",
								   barcode = barcodes)		
		}
		# download the data
		GDCdownload(query_met.hg38)
		# make a data object
		data.hg38 <- GDCprepare_new(query_met.hg38)
		##### revised GDCprepare
		
		# remove downloaded files
		if(removeFiles)unlink(paste0('GDCdata/',"TCGA-",cancer), recursive = TRUE)

		obj = paste0(cancer, '_Meth450')
		assign( obj, data.hg38, envir = .GlobalEnv )
		save(list = obj,file=paste0(folder,obj,'.rda'))
#		return(obj)
		}, error = function(ex) {
  			print(ex);
	})
}

# revised GDCprepare
GDCprepare_new <- function(query,save = FALSE,directory = "GDCdata",summarizedExperiment = TRUE,
               remove.files.prepared = FALSE){
  isServeOK()
  if(missing(query)) stop("Please set query parameter")
  
  test.duplicated.cases <- (any(duplicated(query$results[[1]]$cases)) &
                              !(query$data.type %in% c("Clinical data",
                                                       "Protein expression quantification",
                                                       "Raw intensities",
                                                       "Clinical Supplement",
                                                       "Biospecimen Supplement")))
  
  if(test.duplicated.cases) {
    dup <- query$results[[1]]$cases[duplicated(query$results[[1]]$cases)]
    cols <- c("tags","cases","experimental_strategy","analysis_workflow_type")
    cols <- cols[cols %in% colnames(query$results[[1]])]
    dup <- query$results[[1]][query$results[[1]]$cases %in% dup,cols]
    dup <- dup[order(dup$cases),]
    print(knitr::kable(dup))
    stop("There are samples duplicated. We will not be able to prepare it")
  }
  
  # We save the files in project/source/data.category/data.type/file_id/file_name
  source <- ifelse(query$legacy,"legacy","harmonized")
  files <- file.path(
    query$results[[1]]$project, source,
    gsub(" ","_",query$results[[1]]$data_category),
    gsub(" ","_",query$results[[1]]$data_type),
    gsub(" ","_",query$results[[1]]$file_id),
    gsub(" ","_",query$results[[1]]$file_name)
  )
  
  files <- file.path(directory, files)
  cases <- ifelse(grepl("TCGA|TARGET",query$results[[1]]$project %>% unlist()),query$results[[1]]$cases,query$results[[1]]$sample.submitter_id)
  pb <- txtProgressBar(min = 0, max = length(files), style = 3)
  for (i in seq_along(files)) {
    data <- fread(
      files[i], header = FALSE,
      sep = "\t",
      skip = 1,
      colClasses = c(
        "character", # Composite Element REF
        "numeric" # # beta value
      )
    )
    setnames(data,1,"Composite.Element.REF")
    if (!is.null(cases)) setnames(data,2,cases[i])
    if (i == 1) {
      df <- data
    } else {
      df <- merge(df, data, by = "Composite.Element.REF")
    }
    setTxtProgressBar(pb, i)
  }
#  betas = df
  setDF(df)
  rownames(df) <- df$Composite.Element.REF
  df$Composite.Element.REF <- NULL
  
  # annotate probes with Granges
  colData <-  colDataPrepareTCGA(colnames(df))
  annotation = read.table(file = here('data/HM450.hg38.manifest.tsv'), sep = '\t', header = TRUE)
  rowRanges = annotation[annotation$probeID %in% rownames(df), c("probeID","CpG_chrm","CpG_beg","CpG_end")]
  rownames(rowRanges) = rowRanges$probeID
  rowRanges$probeID = NULL
  colnames(rowRanges)=c("chr","start","end")
  rowRanges = makeGRangesFromDataFrame(rowRanges[complete.cases(rowRanges) & !(rowRanges$chr %in% c("chrX","chrY","chrM")),])
  assay = data.matrix(df[rownames(df) %in% names(rowRanges),])
  colnames(assay) <- rownames(colData)
  rownames(assay) <- as.character(names(rowRanges))
  se = SummarizedExperiment(assays=assay,
                                  rowRanges=rowRanges,
                                     colData=colData)
  return(se)
}


colDataPrepareTCGA <- function(barcode){
  # For the moment this will work only for TCGA Data
  # We should search what TARGET data means
  
  code <- c('01','02','03','04','05','06','07','08','09','10','11',
            '12','13','14','20','40','50','60','61')
  shortLetterCode <- c("TP","TR","TB","TRBM","TAP","TM","TAM","THOC",
                       "TBM","NB","NT","NBC","NEBV","NBM","CELLC","TRB",
                       "CELL","XP","XCL")
  
  definition <- c("Primary solid Tumor", # 01
                  "Recurrent Solid Tumor", # 02
                  "Primary Blood Derived Cancer - Peripheral Blood", # 03
                  "Recurrent Blood Derived Cancer - Bone Marrow", # 04
                  "Additional - New Primary", # 05
                  "Metastatic", # 06
                  "Additional Metastatic", # 07
                  "Human Tumor Original Cells", # 08
                  "Primary Blood Derived Cancer - Bone Marrow", # 09
                  "Blood Derived Normal", # 10
                  "Solid Tissue Normal",  # 11
                  "Buccal Cell Normal",   # 12
                  "EBV Immortalized Normal", # 13
                  "Bone Marrow Normal", # 14
                  "Control Analyte", # 20
                  "Recurrent Blood Derived Cancer - Peripheral Blood", # 40
                  "Cell Lines", # 50
                  "Primary Xenograft Tissue", # 60
                  "Cell Line Derived Xenograft Tissue") # 61
  aux <- DataFrame(code = code,shortLetterCode,definition)
  
  # in case multiple equal barcode
  regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                  "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
  samples <- str_match(barcode,regex)[,1]
  
  ret <- DataFrame(barcode = barcode,
                   patient = substr(barcode, 1, 12),
                   sample = substr(barcode, 1, 16),
                   code = substr(barcode, 14, 15))
  ret <- merge(ret,aux, by = "code", sort = FALSE)
  ret <- ret[match(barcode,ret$barcode),]
  rownames(ret) <- gsub("\\.","-",make.names(ret$barcode,unique=TRUE))
  ret$code <- NULL
  return(DataFrame(ret))
}

# removes specified probes
removeProbes = function(bet, probesToRemove = NULL)
{
	probeWOchrXY=setdiff(rownames(bet),probesToRemove)
	# remove X and Y probes
	dat = bet[probeWOchrXY,]
	# remove all NA probes
	p1 = apply(is.na(dat),1,all)
	noNAprobes = names(p1)[!p1]
	dat = dat[noNAprobes,]
	# remove all non cg probes
	ind = grep('cg.',rownames(dat))
	dat = dat[ind,]

#	colnames(dat) = substr(colnames(dat),1,15)
	return(dat)
}

# returns a list of probes to remove (XY probes, SNPs)
getProbesToExclude = function()
{
	library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
	# get probe locations
	loc = data.frame(Locations)
	# select probes on X and Y chromasomes
	probeChrXY = rownames(loc)[which(loc[,'chr'] %in% c('chrX', 'chrY'))]
	# get annotation of probes to SNPs
	snps = data.frame(SNPs.Illumina)
	snpProbes = rownames(snps)[which(snps[,1] != '' | snps[,1] != '')]
	return(union(probeChrXY,snpProbes))
}

createPhenotype = function(sampleIDs)
{
	n = length(sampleIDs);
	type = vector(length=n);
	tissue = vector(length=n);
	for (i in 1:n)
	{
		t = substr(sampleIDs[i],14,15);
		tissue[i] = t
		# tumor sample
		if (t < 10)	type[i] = 'T';
		if (t>=10 && t<20) type[i] = 'N';
		if (t==20) type[i] = 'CL';
	}
	result = cbind(type,tissue);
	rownames(result)=sampleIDs;
	return (result)
}
# take a path to *_Meth450 object, removes specified probes, saves updated object and normals
removeProbesAndSaveNorm = function(filePath, normFolder = '', outputFolder = '', ...)
{
	require(tools)
	tryCatch({
		# load object with methylation data
		load(filePath)
		# get name of this object
		obj = file_path_sans_ext(basename(filePath))
		# extract cancer type
		cancer = sapply(strsplit(obj,'_',fixed = T),geti,1)
		print(cancer)
		# get object
		methObj = get(obj)	
		# remove probes 
		betaNoXY = removeProbes(assay(methObj), ...)
		# create phenotype
		pheno = createPhenotype(colnames(betaNoXY))
		tumSamp = rownames(pheno)[pheno[,1]=="T"] # 
		normSamp = rownames(pheno)[pheno[,1]=="N"] #
		print(table(pheno[,1]))
		# create a new object with new set of probes
		objNoXY = SummarizedExperiment(assays=list(betaNoXY[,c(normSamp,tumSamp)]),
                     rowRanges=rowRanges(methObj)[rownames(betaNoXY)], colData=colData(methObj)[c(normSamp,tumSamp),])
		# save updated methylation object
		objNoXYname = paste(cancer, 'MethNoXY',sep='_')
		assign( objNoXYname, objNoXY)
		save(list = objNoXYname,file=paste0(outputFolder, objNoXYname,'.rda'))
		
		if(length(normSamp)>0)
		{
		# create an object with normals
			objNorm = SummarizedExperiment(assays=list(betaNoXY[,normSamp]),
						 rowRanges=rowRanges(methObj)[rownames(betaNoXY)], colData=colData(methObj)[normSamp,])
			# save normals
			objNormName = paste(cancer, 'NormSampNoXY',sep='_')
			assign( objNormName, objNorm)
			save(list = objNormName,file=paste0(normFolder,objNormName,'.rda'))		
		}
	}, error = function(ex) {
  			print(ex);
	})
}
removeNAprobes = function(filePath, outputFolder = '')
{
	require(tools)
	library(SummarizedExperiment)
	tryCatch({
		# load object with methylation data
		load(filePath)
		# get name of this object
		obj = file_path_sans_ext(basename(filePath))
		# get object
		methObj = get(obj)
		# remove probes with at least one NA
		methObj = subset(methObj,subset = (rowSums(is.na(assay(methObj))) == 0))
		
		objNoXYname = paste0(obj,'_noNA')
		assign(objNoXYname, methObj)
		save(list = objNoXYname,file=paste0(outputFolder, objNoXYname,'.rda'))
		
	}, error = function(ex) {
  			print(ex);
	})
}
