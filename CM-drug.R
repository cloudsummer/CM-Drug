####Code of CM-drug

#Function definition
.Path001 = function(eachThreadDir, spid, type, idx) {
  tempInfo = sprintf('%s%s.%s.', eachThreadDir, spid, type)
  return(paste0(tempInfo, idx, '.trds'))
}

silenceRemove = function(paths, check = FALSE) {
  stopifnot(is.character(paths))
  stopifnot(is.logical(check))
  
  if (check) {
    existsVector = file.exists(paths) | dir.exists(paths)
    invisible(suppressWarnings(file.remove(paths[existsVector])))
  } else {
    invisible(suppressWarnings(file.remove(paths)))
  }
}

silenceDoGC = function() {
  # silence garbage collection (gc)
  invisible(gc())
}

SysGetSPID = function() {
  # Get the process
  NODE = unlist(strsplit(Sys.info()['nodename'], '\\.'))[1]
  PID = Sys.getpid()
  tempString = sprintf('%s-%s', NODE, PID)
  return(tempString)
}


#temporary directory

SysGetTempDir = function(useRAM = TRUE, usecache = TRUE) {
  
  currentNode = unlist(strsplit(Sys.info()['nodename'], '\\.'), use.names = FALSE)[1]
  
  
  if (usecache) {
    tempDirPath = '/'
    return(tempDirPath)
  }
  
  if (useRAM) {
    tempDirPath = '/the_temp_path/'
    return(tempDirPath)
    
  } else {

    tempDirPath = '/'
    return(tempDirPath)
  }
}

XZFread = function(file, sep = '\t', header = TRUE) {
  stopifnot(is.character(file))
  require(data.table)
  
  tempOutDir = sprintf('%sXZFread_%s/', SysGetTempDir(useRAM = TRUE, usecache = TRUE), SysGetSPID())
  dir.create(tempOutDir, showWarnings = FALSE, recursive = TRUE)
  
  secondStamp = format(Sys.time(), format = '%s')
  fileID = .FileNameScramble(file)
  tempFilePath = sprintf('%s%s_%s.TMP', tempOutDir, secondStamp, fileID)
  xzCommand = sprintf('xz -d -c %s > %s', file, tempFilePath)
  system(command = xzCommand, ignore.stdout = FALSE, ignore.stderr = TRUE)
  
  #read the file by fast way
  tempDT = fread(file = tempFilePath, sep = sep, header = header, showProgress = FALSE)
  
  #Clean the temp
  file.remove(tempFilePath)
  
  return(tempDT)
}

.FileNameScramble = function(charVec) {
  require(stringi)
  
  charVec = stri_trans_tolower(charVec)
  
  charVec = tail(unlist(strsplit(charVec, split = '/'), use.names = FALSE), n = 1)
  
  charVec = unlist(strsplit(charVec, split = ''), use.names = FALSE)
  
  charVec = paste0(na.omit(match(charVec, letters)), collapse = '')
  
  return(charVec)
}


XZSaveRDS = function(obj, file, threads = 32, compression = 6) {

  stopifnot(is.character(file))
  
  silenceRemove(paths = file, check = TRUE)
  
  xzCommand = sprintf('xz -z -T %s -%s > %s', threads, compression, file)
  
  xzConnection = pipe(description = xzCommand, open = 'wb')
  
  saveRDS(object = obj, file = xzConnection)
  
  close(xzConnection)
}


#Core & Minor gene sets
setwd("/home/xiay/1_home_xiay/Project/PD1/R/42_CM-Drug_test/3_Nat_immuno_reviewer.use/work_directory_1/")

list_CM_gene_set=list()

list_CM_gene_set[[1]] <- 
  c("CD8A","HLA-DRA","HLA-DPA1","HLA-DQB1","CD74","IFNG","HLA-DRB5","HLA-DPB1","HLA-DQA1","HLA-DMA","HLA-DRB1","KLRC1","HLA-DQA2",
    "CTSS","KLRD1","HLA-DOA","CIITA","CD8B","KIR2DL4","KLRC3","CD1B","HLA-DMB","CTSB","CD209","CD4","HLA-B",   
    "MICB","HLA-C","B2M","HLA-F","HLA-E","HLA-A","KLRC2")
names(list_CM_gene_set)[1] <- c("Antigen_Processing_and_Presentation")

list_CM_gene_set[[2]] <- 
  c("SH2D1A","IFNG","GZMB","LCP2","PTPN6","PIK3CD","ITGB2","LCK","FASLG","KLRD1","PRF1","HLA-E","CD247","CD48",   
    "CD244","KIR2DL4","KLRC3","IFNB1","KLRC1","MICB","SH2D1B","FCER1G","HCST","HLA-B","HLA-C","HLA-A","TYROBP","KLRC2")
names(list_CM_gene_set)[2] <- c("NaturalKiller_Cell_Cytotoxicity")

list_CM_gene_set[[3]] <- 
  c("K","CD3G","CD3D","ICOS","GRAP2","IFNG","LCP2","CD3E","PTPN6","PIK3CD","NFKBIA","CARD11","CD247","ITK","PTPRC","PDCD1",
    "CTLA4","PIK3R5","IL2","CD8B","CD8A","CD28","CD4","NFKBIE")
names(list_CM_gene_set)[3] <- c("TCR_Signaling_Pathway")

list_CM_gene_set[[4]] <- 
  c("SIRPG","GPR171","CRTAM","GZMA","LAG3","CTSW","PRF1","NKG7","CCR5","C1QB","GZMB","GZMH","LY9","CD7","LAX1","IL7R",
    "ITK","IL2RB","LCP2","KLRG1","SELL","CD8B","CD8A","GNLY")
names(list_CM_gene_set)[4] <- c("Cytotoxiclty_of_ImmuCellAI")

list_CM_gene_set[[5]] <- unique(c(list_CM_gene_set[[1]],list_CM_gene_set[[2]],list_CM_gene_set[[3]],list_CM_gene_set[[4]]))

names(list_CM_gene_set)[5] <- c("Core_gene")

list_CM_gene_set[[6]] <- 
  c("CXCL11","APOBEC3A","CXCL9","IL2","ORM2","ISG20","TCHHL1","OASL","MX1","IL27","TNF","OAS1","CCL8",
    "ISG15","BST2","MX2","CCL1","TLR7","MUC4","ORM1","DEFA3","CCR3","CHIT1","REG3G","C8G","CD40LG",
    "CCR7","MPO","IDO1","IL22","CXCL10","CCL5","PAEP","GNLY","CD8A","CXCL13","CCL4","CCL7","IFNG",
    "PDCD1","CCL4L2","CCL3L3","CCL3","FASLG","FGR","APOBEC3H","HLA-B","MMP12","TLR8","APOBEC3G","IRF7","SLC29A3",
    "HAMP","CD40","NOD2","CXCL2","CTSS","B2M","IRF5","IL10","MARCO","BPI","HCK","CYBB","IL6", 
    "CXCR1","DEFA4","DEFA1","S100A12","CAMP","PGLYRP1","CCL13","AZU1")
names(list_CM_gene_set)[6] <- c("Antimicrobials")

list_CM_gene_set[[7]] <- 
  c("CD79A","CD19","CD79B","CARD11","CD72","INPP5D","PLCG2","PIK3CG","PTPN6","PRKCB","RAC2","VAV1","BTK","SYK","PIK3CD","PIK3R5",
    "NFKBIA","LYN","CR2","LILRB3","NFKBIE","JUN")
names(list_CM_gene_set)[7] <- c("BCR_Signaling_Pathway")

list_CM_gene_set[[8]] <- unique(c(list_CM_gene_set[[6]],list_CM_gene_set[[7]]))

names(list_CM_gene_set)[8] <- c("Minor_gene")

#we set the name of CM gene sets as the "super" initially. So in the code, for consistency, we keep the name 
super.third.human.pd1.all <- list_CM_gene_set

saveRDS(super.third.human.pd1.all,"./super.third.human.pd1.all.RDS")

################################################################################
################################################################################

require(data.table)
require(stringi)
require(parallel)

# Declaring Global Variables
inDir = c('Data'='Data1')

options(
  stringsAsFactors = FALSE,
  warn = 1
)

# Global setting
# Output
outDir = './Data/Compound_Data/'
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Custom Functions
SelectionResolver = function(conditionDF, controlDF) {
  # Output Resolving Function
  finalSampleSize = min(conditionDF$sampleCount, controlDF$sampleCount, 20)
  tempVector = c(
    'source_dataset' = conditionDF$source_dataset,
    'rna_centre' = conditionDF$rna_centre,
    'cell_id' = conditionDF$cell_id,
    
    'trt_type' = conditionDF$pert_type,
    'trt_cval' = conditionDF$pert_cval,
    'trt_iname' = conditionDF$pert_iname,
    'trt_cdose' = conditionDF$pert_cdose,
    'trt_ctime' = conditionDF$pert_ctime,
    'trt_orig_sampleCount' = conditionDF$sampleCount,
    'trt_final_sampleCount' = finalSampleSize,
    
    'ctl_type' = controlDF$pert_type,
    'ctl_cval' = controlDF$pert_cval,
    'ctl_iname' = controlDF$pert_iname,
    'ctl_cdose' = controlDF$pert_cdose,
    'ctl_ctime' = controlDF$pert_ctime,
    'ctl_orig_sampleCount' = controlDF$sampleCount,
    'ctl_final_sampleCount' = finalSampleSize
  )
  return(tempVector)
}

#Matched the Sample and Control
finalSampleMetadata = NULL
finalMatchedMetadata = NULL

for (currentIndex in 1) {
  currentCase = names(inDir)[currentIndex]
  currentDir = inDir[currentIndex]
  
  cat(sprintf('\nLoading Metadata - %s.\n', currentCase))
  tempPath = sprintf('pertInfo.%s', currentDir)
  metadata = fread(input = tempPath, sep = '\t', header = TRUE, data.table = FALSE)
  
  # ----- Sample Metadata -----
  
  cat('Generating Sample Metadata.\n')
  
  #in original data, -666 means NA, now tidy the data 
  metadata = metadata[!grepl('^-666', metadata$pert_iname), ]
  metadata$pert_dose_unit[metadata$pert_dose == -666] = NA
  metadata$pert_dose[metadata$pert_dose == -666] = NA
  metadata$pert_dose[metadata$pert_time == -666] = NA
  
  metadata$pert_cdose = sprintf(
    '%s%s',
    metadata$pert_dose,
    metadata$pert_dose_unit
  )
  metadata$pert_cdose[grepl('^NA', metadata$pert_cdose)] = ''
  
  metadata$pert_ctime = sprintf(
    '%s%s',
    metadata$pert_time,
    metadata$pert_time_unit
  )
  
  metadata$pert_cval = sprintf(
    '%s %s for %s in %s',
    metadata$pert_cdose,
    metadata$pert_iname,
    metadata$pert_ctime,
    metadata$cell_id
  )
  metadata$pert_cval = stri_trim_left(metadata$pert_cval)
  
  metadata$rna_centre = sapply(strsplit(x = as.character(metadata$rna_plate), split = '_'), 
                               USE.NAMES = FALSE, FUN = function(i) {return(i[[1]])})
  
  selectionVector = c(
    'Num',
    'rna_plate',
    'rna_centre',
    'pert_iname',
    'pert_type',
    'pert_dose',
    'pert_dose_unit',
    'pert_time',
    'pert_cdose',
    'pert_ctime',
    'pert_cval',
    'cell_id'
  )
  sampleMetadata = metadata[, selectionVector]
  sampleMetadata$source_dataset = rep(currentCase, nrow(metadata))
  
  # Variable Cleanup
  rm(selectionVector, metadata)
  
  # ----- Generate Baseline-Matched Condition Metadata -----
  
  cat('Generating Baseline-Matched Condition Metadata.\n')
  
  # Match Contrast-Baseline
  selectionVector = c(
    'source_dataset',
    'pert_cval',
    'rna_centre',
    'pert_type',
    'cell_id',
    'pert_iname',
    'pert_ctime',
    'pert_cdose'
  )
  focusDF = unique(sampleMetadata[, selectionVector])
  
  #Sample Counts
  sampleMetadata = as.data.table(sampleMetadata)
  setkey(sampleMetadata, rna_centre, pert_type, pert_cval)
  focusDF$sampleCount = unlist(mclapply(1:nrow(focusDF), mc.preschedule = TRUE, mc.cores = 64, mc.cleanup = TRUE, FUN = function(i) {
    tempCondition = focusDF[i, ]
    return(sampleMetadata[.(tempCondition$rna_centre, tempCondition$pert_type, tempCondition$pert_cval), .N, nomatch = 0])
  }), recursive = FALSE, use.names = FALSE)
  sampleMetadata = as.data.frame(sampleMetadata)
  
  # Filter out Sample Size == 1
  focusDF = subset(focusDF, sampleCount > 1)
  
  # Prepare Subsets
  casesDF = subset(focusDF, grepl('^trt', focusDF$pert_type))
  controlsDF = subset(focusDF, grepl('^ctl', focusDF$pert_type))
  
  tempList = mclapply(1:nrow(casesDF), mc.preschedule = TRUE, mc.cores = 32, mc.cleanup = TRUE, FUN = function(currentRowIndex) {
    currentCondition = casesDF[currentRowIndex, ]
    currentCentre = currentCondition$rna_centre
    currentType = currentCondition$pert_type
    currentCell = currentCondition$cell_id
    currentTime = currentCondition$pert_ctime
    currentDose = currentCondition$pert_cdose
    
    controlSubset = subset(controlsDF, rna_centre == currentCentre & cell_id == currentCell & pert_ctime == currentTime)
    
    # Compound Data Resolving
    if (currentType == 'trt_cp') {
      controlSubset = subset(controlSubset, pert_type %in% c('ctl_vehicle', 'ctl_untrt'))
      
      # SKIP if No Matching: 
      if (nrow(controlSubset) == 0) {
        return(NULL)
      }
      
      if ('DMSO' %in% controlSubset$pert_iname) {
        # Baseline: DMSO
        controlSubset = subset(controlSubset, pert_iname == 'DMSO')
        finalOutput = controlSubset[order(controlSubset$sampleCount, decreasing = TRUE), ][1, ]
        return(SelectionResolver(conditionDF = currentCondition, controlDF = finalOutput))
      }
      
      if (any(c('PBS', 'H2O', 'UnTrt') %in% controlSubset$pert_iname)) {
        # Baseline: PBS, H2O, Untreated Cells
        controlSubset = subset(controlSubset, pert_iname %in% c('PBS', 'H2O', 'UnTrt'))
        finalOutput = controlSubset[order(controlSubset$sampleCount, decreasing = TRUE), ][1, ]
        return(SelectionResolver(conditionDF = currentCondition, controlDF = finalOutput))
        
      } else {
        # Baseline of Last Resort
        finalOutput = controlSubset[order(controlSubset$sampleCount, decreasing = TRUE), ][1, ]
        return(SelectionResolver(conditionDF = currentCondition, controlDF = finalOutput))
      }
    }
    
  })
  
  # Assemble Condition Metadata
  cat('Assemble Condition Metadata.\n')
  matchedMetadata = as.data.frame(do.call('rbind', tempList))
  
  # Append Metadata to Finalized Version
  finalSampleMetadata = rbind(finalSampleMetadata, sampleMetadata)
  finalMatchedMetadata = rbind(finalMatchedMetadata, matchedMetadata)
  
  # Variable Cleanup
  rm(currentCase, currentDir)
  rm(selectionVector, focusDF, casesDF, controlsDF, tempList)
  rm(sampleMetadata, matchedMetadata)
}

# Appending Condition ID
finalMatchedMetadata = data.frame(
  case_ID = paste0('CM.', 1:nrow(finalMatchedMetadata)),
  finalMatchedMetadata
)

#Switching type
finalMatchedMetadata$trt_orig_sampleCount = as.numeric(finalMatchedMetadata$trt_orig_sampleCount)
finalMatchedMetadata$trt_final_sampleCount = as.numeric(finalMatchedMetadata$trt_final_sampleCount)

finalMatchedMetadata$ctl_orig_sampleCount = as.numeric(finalMatchedMetadata$ctl_orig_sampleCount)
finalMatchedMetadata$ctl_final_sampleCount = as.numeric(finalMatchedMetadata$ctl_final_sampleCount)

#Saving
cat('Now saving...\n')

#Sample Metadata
tempPath = sprintf('%ssample.metadata.RDS.XZ', outDir)
XZSaveRDS(obj = finalSampleMetadata, file = tempPath)

#Matched Metadata
tempPath = sprintf('%scondition.metadata.RDS.XZ', outDir)
XZSaveRDS(obj = finalMatchedMetadata, file = tempPath)

# cat(sprintf('START TIME: %s\n', .startTime))
# START TIME: Fri Apr  1 10:08:28 2022
# cat(sprintf('END TIME: %s\n\n', date()))
# END TIME: Fri Apr  1 20:23:49 2022

# Load Libraries
require(cmapR)
require(data.table)
require(foreach)
require(doParallel)

# Declaring Global Variables
inPath = c('Data.1'='Data_A549.1.gctx'
)
metadataDir = './Data/Compound_Data/'

options(
  stringsAsFactors = FALSE,
  warn = 1
)

#Preparation 
################################################################################
# Creation of Output Directory
outDir = 'Data/Compound_Data/'
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Creation of Interim RAM-backed Output Directory
ramDir = '/home/xiay/1_home_xiay/ramdisk/'
dir.create(ramDir, recursive = TRUE, showWarnings = FALSE)

################################################################################
cat('Loading Metadata.\n')
tempPath = sprintf('%ssample.metadata.RDS.XZ', metadataDir)
sampleDF = readRDS(tempPath)

tempPath = sprintf('%scondition.metadata.RDS.XZ', metadataDir)
matchedDF = readRDS(tempPath)

# Sample Metadata Trimming
uniqueCVals = unique(c(matchedDF$trt_cval, matchedDF$ctl_cval))
sampleDF = subset(sampleDF, pert_cval %in% uniqueCVals)
rm(uniqueCVals)

# Enable data.table optimization
sampleDF = as.data.table(sampleDF)
setkey(sampleDF, rna_centre, pert_type, pert_cval)

cat('Loading Data.\n')
referenceRow = cmapR::read_gctx_ids(gctx_path = inPath['Data.1'],
                                    #dimension = 'row'
)
dataList = lapply(inPath, function(tempPath) {
  tempMatrix = parse_gctx(fname = tempPath)@mat[referenceRow, ]
  silenceDoGC()
  return(tempMatrix)
})

gctxMatrix = dataList[[1]]
rm(dataList)
silenceDoGC()

# Data Matrix Trimming
gctxMatrix = gctxMatrix[, sampleDF$Num]
silenceDoGC()

# Calculating
cat('Calculating...\n')

# Initiate Multi-Thread
.startTime = date()
currentSPID = SysGetSPID()
referenceColumn = matchedDF$case_ID
caseCount = length(matchedDF$case_ID)
geneCount = length(referenceRow)
registerDoParallel(cores = 180)

# Cleaning Cache
allPath = c(
  .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_fc', idx = 1:caseCount)
)
silenceRemove(paths = allPath, check = TRUE)

# Parallel Processing
cat('Begin Calculating Job: ')
tempOutput = foreach(i = 1:caseCount, .inorder = TRUE) %dopar% {
  # Status Update
  if (i %% 500 == 0) {
    cat(i, ' . ', sep = '')
  }
  
  # Obtain Matching-Sample Mapping
  currentTask = matchedDF[i, ]
  trtSamples = sampleDF[.(currentTask$rna_centre, currentTask$trt_type, currentTask$trt_cval), Num]
  ctlSamples = sampleDF[.(currentTask$rna_centre, currentTask$ctl_type, currentTask$ctl_cval), Num]
  
  # Sample-Size Control
  if (currentTask$trt_orig_sampleCount != currentTask$trt_final_sampleCount) {
    tempLogical = sample(x = 1:currentTask$trt_orig_sampleCount, size = currentTask$trt_orig_sampleCount, replace = FALSE)
    tempLogical = (tempLogical %in% 1:currentTask$trt_final_sampleCount)
    trtSamples = trtSamples[tempLogical]
    rm(tempLogical)
  }
  
  if (currentTask$ctl_orig_sampleCount != currentTask$ctl_final_sampleCount) {
    tempLogical = sample(x = 1:currentTask$ctl_orig_sampleCount, size = currentTask$ctl_orig_sampleCount, replace = FALSE)
    tempLogical = (tempLogical %in% 1:currentTask$ctl_final_sampleCount)
    ctlSamples = ctlSamples[tempLogical]
    rm(tempLogical)
  }
  
  trtData = gctxMatrix[, trtSamples]
  ctlData = gctxMatrix[, ctlSamples]
  
  tempFC = rowMeans(trtData) - rowMeans(ctlData)

  # Save File
  tempPath = .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_fc', idx = i)
  XZSaveRDS(obj = tempFC, file = tempPath)
  
  # Variable Cleanup
  rm(currentTask, trtSamples, ctlSamples)
  rm(trtData, ctlData)
  rm(tempResult, tempFC)
  silenceDoGC()
  return(NULL)
}

# stop the multi-threads
registerDoSEQ()

#clean-up the variable 
rm(matchedDF, sampleDF, gctxMatrix, allPath, tempOutput)
silenceDoGC()
cat('clean-up have been done.\n')

tempMatrix = foreach(i = 1:caseCount, .inorder = TRUE, .combine = cbind, .maxcombine = 1000) %do% {
  tempPath = .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_fc', idx = i)
  return(readRDS(tempPath))
}
colnames(tempMatrix) = referenceColumn
rownames(tempMatrix) = referenceRow

tempPath = sprintf('%sfc.matrix.RDS.XZ', outDir)
XZSaveRDS(obj = tempMatrix, file = tempPath)
rm(tempMatrix)
silenceDoGC()
silenceRemove(paths = .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_fc', idx = 1:caseCount))

# Directory Cleanup
silenceRemove(ramDir)

# Print Timestamp
cat(sprintf('START TIME: %s\n', .startTime))
cat(sprintf('END TIME: %s\n\n', date()))

# cat(sprintf('START TIME: %s\n', .startTime))
# START TIME: Fri Apr  1 08:08:28 2022
# cat(sprintf('END TIME: %s\n\n', date()))
# END TIME: Fri Apr  1 23:23:49 2022
# stopCluster(cl)
################################################################################
################################################################################



fc.matrix <- readRDS("./Data/Compound_Data/fc.matrix.RDS.XZ")
fc.compound_pert.gtc<- new("GCT", mat=fc.matrix)
write_gctx(fc.compound_pert.gtc, 
           #compression_level = 9,
           "./Data/Compound_Data/fc.compound_pert")

rm(list = ls())
setwd("/home/xiay/1_home_xiay/Project/PD1/R/42_CM-Drug_test/3_Nat_immuno_reviewer.use/work_directory_1/")
metadata.condition <- readRDS("./Data/Compound_Data/condition.metadata.RDS.XZ")


data_trt_cp <- list()
data_trt_cp.df <- list()
fgsea.sam.trt_cp <- list()
fgsea.res.trt_cp <- list()
trt_cp_number.group <- list()

trt_cp_number <- which(metadata.condition$trt_type == "trt_cp")

trt_cp_number.group[[1]] <- trt_cp_number[1:length(trt_cp_number)]

template <- parse_gctx("Data/Compound_Data/fc.compound_pert_n20936x12328.gctx",
                       rid=1:12328, cid=1:10)@mat %>% as.data.frame()

template <- template %>% dplyr::mutate(gene_id=rownames(template ))

gene_df <- read.delim("./gene_df.1.txt")

template$gene_id <- as.integer(template$gene_id)

template <- dplyr::left_join(template,gene_df,by= "gene_id")

saveRDS(template,"template.RDS")
save.image("template.RData")

################################################################################
################################################################################
#Use fgsea to perform GSEA 
options(
  stringsAsFactors = FALSE,
  warn = 0
)
library(fgsea)
library(cmapR)
library(tidyverse)
setwd("~/1_home_xiay/Project/PD1/R/42_CM-Drug_test/3_Nat_immuno_reviewer.use/work_directory_1/")
resultDir="./result/cp/"
dir.create(resultDir, recursive = TRUE, showWarnings = FALSE)
load("./template.RData")

#library(tictoc)
library(furrr)
library(future)
plan(multisession, workers = 120)

.startTime = date()

for(j in c(1)){
  setwd("~/1_home_xiay/Project/PD1/R/42_CM-Drug_test/3_Nat_immuno_reviewer.use/work_directory_1/")
  load("./template.RData")
  
  fun_cmap_fgsea <- function(x){
    x <- x %>% unlist() 
    names(x) <-template$gene_symbol
    fgsea.res <- fgsea(pathways = genesets, stats = x,eps= 0.0, minSize  = 5, maxSize  = 500)
    return(fgsea.res)
  }
  genesets <- readRDS("super.third.human.pd1.all.RDS")
  data_trt_cp[[j]] <- parse_gctx("./Data/Compound_Data/fc.compound_pert_n20936x12328.gctx", 
                                 rid=1:12328, cid=trt_cp_number.group[[j]])
  data_trt_cp.df[[j]] <- as.data.frame(data_trt_cp[[j]]@mat)
  fgsea.res.trt_cp[[j]] <- furrr::future_map(data_trt_cp.df[[j]], ~ fun_cmap_fgsea(.x))
  setwd("./result/cp/");saveRDS(fgsea.res.trt_cp[[j]],paste("c",j,"_fgsea.res.trt_cp_super.RDS",sep = ""));rm(list=ls());
}

cat(sprintf('START TIME: %s\n', .startTime))
cat(sprintf('END TIME: %s\n\n', date()))


################################################################################
################################################################################
#tidy the results

setwd("/home/xiay/1_home_xiay/Project/PD1/R/42_CM-Drug_test/3_Nat_immuno_reviewer.use/work_directory_1/result/cp/")
load("/home/xiay/1_home_xiay/Project/PD1/R/42_CM-Drug_test/3_Nat_immuno_reviewer.use/work_directory_1/template.RData")

fgsea.res.trt_cp <- list()
for(i in 1){fgsea.res.trt_cp[[i]] <- readRDS(paste("./c",i,"_fgsea.res.trt_cp_super.RDS",sep = ""))}

fgsea.res.tidy <- list()

for(i in 1:8){print(i);
  fun_tidy= function(GSEA_res_list){
    x <- furrr::future_map_dfr(GSEA_res_list, ~ .x[i,])
    return(x)
  }
  fgsea.res.tidy[[i]] <- furrr::future_map(fgsea.res.trt_cp,~ fun_tidy(.x)) %>% Reduce(rbind,.)
}

fgsea.clue.order_with_id <- list()
for(i in 1:8){fgsea.clue.order_with_id[[i]] <- fgsea.res.tidy[[i]] %>% dplyr::mutate(id=rownames(fgsea.res.tidy[[i]]))}

for(i in 1:8){
  names(fgsea.clue.order_with_id)[i] <- fgsea.clue.order_with_id[[i]][1,1]
}

saveRDS(fgsea.clue.order_with_id,"fgsea.clue.order_with_id.RDS")

setwd("/home/xiay/1_home_xiay/Project/PD1/R/42_CM-Drug_test/3_Nat_immuno_reviewer.use/work_directory_1/")

fgsea.clue.order <- readRDS("./result/cp/fgsea.clue.order_with_id.RDS")

#screening strategy1
# fgsea.clue.screen1 <- purrr::map(fgsea.clue.order,~ dplyr::filter(.x,padj<0.05,NES>0))
# fgsea.clue.screen2 <- purrr::map(fgsea.clue.screen1,~  .x$id)
# fgsea.clue.screen3 <- Reduce(intersect,fgsea.clue.screen2)
# good.id <- fgsea.clue.screen3

#screening strategy2
#Considering that certain drugs significantly enhance specific pathways but show slightly insufficient significance in one or two other pathways, to prevent the oversight of such drugs, we have implemented an alternative threshold screening allowing for some flexibility. Specifically, we allow an adjusted p-value between 0.05 and 0.2 for enrichment results in two pathways, while the adjusted p-value for other enriched pathways must be less than 0.05.
fgsea.clue.screen1 <- purrr::map(fgsea.clue.order,~ dplyr::filter(.x,padj<0.2,NES>0))
fgsea.clue.screen2 <- purrr::map(fgsea.clue.screen1,~  .x$id)
fgsea.clue.screen3 <- Reduce(intersect,fgsea.clue.screen2)
fgsea.clue.screen4 <- purrr::map(fgsea.clue.order,~ .x[as.integer(fgsea.clue.screen3),])
fgsea.clue.screen5 <- purrr::map(fgsea.clue.screen4 ,~ .x[,3])
fgsea.clue.screen6 <- Reduce(cbind,fgsea.clue.screen5)
fgsea.clue.screen6 <- cbind(fgsea.clue.screen6,fgsea.clue.screen4[[1]]$id)
fgsea.clue.screen6 <- fgsea.clue.screen6 %>% as.data.frame()
fgsea.clue.screen6.row <- purrr::map(as.data.frame(t(fgsea.clue.screen6[,1:8])), ~ .x)
fgsea.clue.screen6.row <- purrr::map(fgsea.clue.screen6.row,~ as.numeric(.))
fgsea.clue.screen6.count <- purrr::map(fgsea.clue.screen6.row,function(x){count <- 0;
for(i in 1:8){if(x[i]<0.05){count <- count+1}};return(count);count <- 0})
unlist(fgsea.clue.screen6.count) %>% as.data.frame() -> temp
colnames(fgsea.clue.screen6)[1:8] <- paste(colnames(fgsea.clue.screen6)[1],1:8,sep = "");
colnames(fgsea.clue.screen6)[9] <- "id"
fgsea.clue.screen6 <- fgsea.clue.screen6  %>% mutate(num_of_0.05=temp$.)
fgsea.clue.screen.end <- dplyr::filter(fgsea.clue.screen6,num_of_0.05>=7)
fgsea.clue.screen.id <- fgsea.clue.screen.end$id
good.id <- fgsea.clue.screen.id

##screen data
load("./template.RData")

list.good <- list()
for(i in 1:8){list.good[[i]] <- fgsea.clue.order[[i]][as.numeric(good.id),]}
good.df <- as.data.frame(matrix(data = NA,nrow =dim(list.good[[1]])[1],ncol = 8));{for(i in 1:8){good.df[,i] <- list.good[[i]][,6]}}
colnames(good.df) <- names(fgsea.clue.order)

good.score <- vector();
################################################################################
#CM-Score
for(i in 1:dim(list.good[[1]])[1]){good.score[i] <- 0.4986+0.0969*good.df[,1][i]+0.0892*good.df[,2][i]+0.0307*good.df[,3][i]+0.0117*good.df[,4][i]+0.0124*good.df[,6][i]+ 0.0213*good.df[,7][i]}

metadata.condition <- readRDS("./Data/Compound_Data/condition.metadata.RDS.XZ")

good.df.withscore <- data.frame(good.score,good.df)

good.df.withscore <- data.frame(rownames(good.df.withscore),good.df.withscore)

colnames(good.df.withscore)[1] <- "good_id_index"

good.df.withscore.order <- good.df.withscore %>% arrange(desc(good.score))

good.id.order <- good.id[as.numeric(good.df.withscore.order$good_id_index)]

good.id.order.trt_number <- trt_cp_number[as.numeric(good.id.order)]

sig_info.choose <- metadata.condition[good.id.order.trt_number,]

sig_info.choose <- data.frame(1:dim(list.good[[1]])[1],sig_info.choose)

colnames(sig_info.choose)[1] <- "id"

sig_info.choose <- sig_info.choose %>% mutate(Pert_Score=good.df.withscore.order$good.score)

sig_info.choose.fgsea <- cbind(sig_info.choose,good.df.withscore.order)

################################################################################
sig_info.choose.A549<- sig_info.choose.fgsea %>% dplyr::filter(cell_id == "A549")

drug.choose.detail <- table(sig_info.choose$trt_iname) %>% as.data.frame()%>% arrange(desc(Freq))

drug.choose.detail <- drug.choose.detail %>% mutate( pert_score_max=NA,pert_score_detail=NA)

colnames(drug.choose.detail)[1] <- "Compound"

colnames(drug.choose.detail)[2] <- "Times in A549"

drug.choose.detail$pert_score_max <-purrr::map(drug.choose.detail$Compound,~ dplyr::filter(sig_info.choose,trt_iname==.x) 
                                               %>% select(Pert_Score)%>% max() %>% round(2) %>% unlist())

drug.choose.detail$pert_score_max<- drug.choose.detail$pert_score_max %>% unlist()

drug.choose.detail$pert_score_detail <-purrr::map(drug.choose.detail$Compound,~ dplyr::filter(sig_info.choose,trt_iname==.x) %>% 
                                                    select(Pert_Score)%>% round(2) %>% unlist(.) %>% as.data.frame() %>% .[,1])

fun_paste <- function(x,y){
  paste(x,y,sep="_")
}

drug.choose.detail$pert_score_detail <- purrr::map(drug.choose.detail$pert_score_detail,~ Reduce(fun_paste,.))

drug.choose.detail <- purrr::map_df(drug.choose.detail,~ unlist(.)) 

drug.choose.detail[is.na( drug.choose.detail)] <- -666 #in original metadata -666 means NA

#save in excel file-type
library(openxlsx)
openxlsx::write.xlsx(x = drug.choose.detail , file = "./result/results_with_CM_Drug.2.xlsx",
                     sheetName = "screenResult", rownames = FALSE)
