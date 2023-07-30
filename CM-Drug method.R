#### Use LINCS data to screen ICB synergistic boosters by the CM-Drug method
#### The data needed is downloaded form LINCS(CLUE), https://clue.io/releases/data-dashboard

#set path and global function
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
  # Obtain Server-Process Combined ID
  NODE = unlist(strsplit(Sys.info()['nodename'], '\\.'))[1]
  PID = Sys.getpid()
  tempString = sprintf('%s-%s', NODE, PID)
  return(tempString)
}


#temporary directory (for save the usage of RAM and disk)

SysGetTempDir = function(useRAM = TRUE, usecache = TRUE) {
  
  currentNode = unlist(strsplit(Sys.info()['nodename'], '\\.'), use.names = FALSE)[1]
  
  
  if (usecache) {
    tempDirPath = '/'
    return(tempDirPath)
  }
  
  # Check RAM
  if (useRAM) {
    tempDirPath = '/the_temp_path/'
    return(tempDirPath)
    
  } else {
    # Revert to output-direcotry
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
  
  # Fast-Read
  tempDT = fread(file = tempFilePath, sep = sep, header = header, showProgress = FALSE)
  
  # Cleanup
  file.remove(tempFilePath)
  
  return(tempDT)
}

.FileNameScramble = function(charVec) {
  # Hidden Function For Scrambling File Names
  require(stringi)
  
  charVec = stri_trans_tolower(charVec)
  charVec = tail(unlist(strsplit(charVec, split = '/'), use.names = FALSE), n = 1)
  charVec = unlist(strsplit(charVec, split = ''), use.names = FALSE)
  charVec = paste0(na.omit(match(charVec, letters)), collapse = '')
  return(charVec)
}


XZSaveRDS = function(obj, file, threads = 32, compression = 6) {
  # Parallel XZ-Compression SaveRDS
  stopifnot(is.character(file))
  
  silenceRemove(paths = file, check = TRUE)
  
  xzCommand = sprintf('xz -z -T %s -%s > %s', threads, compression, file)
  xzConnection = pipe(description = xzCommand, open = 'wb')
  saveRDS(object = obj, file = xzConnection)
  close(xzConnection)
}


#Core & Minor gene sets

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
names(list_CM_gene_set)[4] <- c("$Cytotoxiclty_of_ImmuCellAI")

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

saveRDS(super.third.human.pd1.all,"../GSEA_Genesets/super.third.human.pd1.all.RDS")
################################################################################
################################################################################

# Load Libraries
require(data.table)
require(stringi)
require(parallel)

#make tidy inst_info of key information

setwd("./1_home_xiay/Project/PD1/R/35_LINCS_2020/")

instinfo_beta <- fread(input = "../../LINCS/3_LINCS_2020/level3/instinfo_beta.txt", sep = '\t', header = TRUE, data.table = FALSE)

Inst_diy <- as.data.frame(matrix(data=NA,nrow=3026460,ncol = 11))
colnames(Inst_diy) <- colnames(Inst_Info)
Inst_diy$inst_id <- instinfo_beta$sample_id
Inst_diy$rna_plate <- instinfo_beta$rna_plate
Inst_diy$rna_well <- instinfo_beta$rna_well
Inst_diy$pert_id <- instinfo_beta$pert_id
Inst_diy$pert_iname<- instinfo_beta$cmap_name
Inst_diy$pert_type<- instinfo_beta$pert_type
Inst_diy$pert_dose<- instinfo_beta$pert_dose
Inst_diy$pert_dose_unit<- instinfo_beta$pert_dose_unit
Inst_diy$pert_time<- instinfo_beta$pert_time
Inst_diy$pert_time_unit<- instinfo_beta$pert_time_unit
Inst_diy$cell_id<- instinfo_beta$cell_iname
#write.table(Inst_diy,"../../LINCS/3_LINCS_2020/level3/Inst_diy.txt",sep = "\t",col.names = TRUE,row.names = FALSE)

# Declaring Global Variables
inDir = c('LINCS_2020'='../../LINCS/3_LINCS_2020/level3/'
)

options(
  stringsAsFactors = FALSE,
  warn = 1
)

# Global setting
# Output
outDir = 'Data/LINCS/Level3/Compound_Data/'
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Custom Functions
SelectionResolver = function(conditionDF, controlDF) {
  # Output Resolving Function
  finalSampleSize = min(conditionDF$sample_size, controlDF$sample_size, 20)
  tempVector = c(
    'source_dataset' = conditionDF$source_dataset,
    'rna_centre' = conditionDF$rna_centre,
    'cell_id' = conditionDF$cell_id,
    
    'trt_type' = conditionDF$pert_type,
    'trt_cval' = conditionDF$pert_cval,
    'trt_iname' = conditionDF$pert_iname,
    'trt_cdose' = conditionDF$pert_cdose,
    'trt_ctime' = conditionDF$pert_ctime,
    'trt_orig_sample_size' = conditionDF$sample_size,
    'trt_final_sample_size' = finalSampleSize,
    
    'ctl_type' = controlDF$pert_type,
    'ctl_cval' = controlDF$pert_cval,
    'ctl_iname' = controlDF$pert_iname,
    'ctl_cdose' = controlDF$pert_cdose,
    'ctl_ctime' = controlDF$pert_ctime,
    'ctl_orig_sample_size' = controlDF$sample_size,
    'ctl_final_sample_size' = finalSampleSize
  )
  return(tempVector)
}

# Iterative Processing -------------------------------------------------------
# Container For Final Metadata (Sample and Baseline-Matched)
finalSampleMetadata = NULL
finalMatchedMetadata = NULL

for (currentIndex in 1) {
  currentCase = names(inDir)[currentIndex]
  currentDir = inDir[currentIndex]
  
  cat(sprintf('\nLoading Metadata - %s.\n', currentCase))
  tempPath = sprintf('%sInst_diy.txt', currentDir)
  metadata = fread(input = tempPath, sep = '\t', header = TRUE, data.table = FALSE)
  
  # ----- Sample Metadata -----
  
  cat('Generating Sample Metadata.\n')
  
  #in LINCS data, -666 means NA, now tidy the data 
  metadata = metadata[!grepl('^-666', metadata$pert_iname), ]
  metadata$pert_dose_unit[metadata$pert_dose == -666] = NA
  metadata$pert_dose[metadata$pert_dose == -666] = NA
  metadata$pert_dose[metadata$pert_time == -666] = NA
  
  #Column Names fixed to make it more readable
  if (currentIndex == 2) {
    colnames(metadata)[colnames(metadata) == 'det_plate'] = 'rna_plate'
    colnames(metadata)[colnames(metadata) == 'det_well'] = 'rna_well'
  }
  
  # Generate Combination Identifiers
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
  
  # Generate "Experiment Centre" Identifiers
  metadata$rna_centre = sapply(strsplit(x = as.character(metadata$rna_plate), split = '_'), USE.NAMES = FALSE, FUN = function(i) {return(i[[1]])})
  
  # Generate Final Sample-Level Metadata
  selectionVector = c(
    'inst_id',
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
  
  # Populating Sample Counts
  sampleMetadata = as.data.table(sampleMetadata)
  setkey(sampleMetadata, rna_centre, pert_type, pert_cval)
  focusDF$sample_size = unlist(mclapply(1:nrow(focusDF), mc.preschedule = TRUE, mc.cores = 64, mc.cleanup = TRUE, FUN = function(i) {
    tempCondition = focusDF[i, ]
    return(sampleMetadata[.(tempCondition$rna_centre, tempCondition$pert_type, tempCondition$pert_cval), .N, nomatch = 0])
  }), recursive = FALSE, use.names = FALSE)
  sampleMetadata = as.data.frame(sampleMetadata)
  
  # Filter out Sample Size == 1
  focusDF = subset(focusDF, sample_size > 1)
  
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
        finalOutput = controlSubset[order(controlSubset$sample_size, decreasing = TRUE), ][1, ]
        return(SelectionResolver(conditionDF = currentCondition, controlDF = finalOutput))
      }
      
      if (any(c('PBS', 'H2O', 'UnTrt') %in% controlSubset$pert_iname)) {
        # Baseline: PBS, H2O, Untreated Cells
        controlSubset = subset(controlSubset, pert_iname %in% c('PBS', 'H2O', 'UnTrt'))
        finalOutput = controlSubset[order(controlSubset$sample_size, decreasing = TRUE), ][1, ]
        return(SelectionResolver(conditionDF = currentCondition, controlDF = finalOutput))
        
      } else {
        # Baseline of Last Resort
        finalOutput = controlSubset[order(controlSubset$sample_size, decreasing = TRUE), ][1, ]
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
  case_ID = paste0('CC.', 1:nrow(finalMatchedMetadata)),
  finalMatchedMetadata
)

#Switching type
finalMatchedMetadata$trt_orig_sample_size = as.numeric(finalMatchedMetadata$trt_orig_sample_size)
finalMatchedMetadata$trt_final_sample_size = as.numeric(finalMatchedMetadata$trt_final_sample_size)

finalMatchedMetadata$ctl_orig_sample_size = as.numeric(finalMatchedMetadata$ctl_orig_sample_size)
finalMatchedMetadata$ctl_final_sample_size = as.numeric(finalMatchedMetadata$ctl_final_sample_size)

#Saving
cat('Saving...\n')

# Metadata of the samples
tempPath = sprintf('%ssample.metadata.RDS.XZ', outDir)# LINCS Level 3: Differential Expression

# Load Libraries
require(cmapR)
require(data.table)
require(foreach)
require(doParallel)

# Declaring Global Variables
inPath = c('LINCS_2020'='../../LINCS/3_LINCS_2020/level3/level3_beta_all_n3026460x12328.gctx'
)
metadataDir = 'Data/LINCS/Level3/Compound_Data/'

options(
  stringsAsFactors = FALSE,
  warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
outDir = 'Data/LINCS/Level3/Compound_Data/'
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Creation of Interim RAM-backed Output Directory
ramDir = '/home/xiay/1_home_xiay/ramdisk/'
dir.create(ramDir, recursive = TRUE, showWarnings = FALSE)

# Main Processing -------------------------------------------------------
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
referenceRow = cmapR::read_gctx_ids(gctx_path = inPath['LINCS_2020'],
                                    #dimension = 'row'
)
dataList = lapply(inPath, function(tempPath) {
  tempMatrix = parse.gctx(fname = tempPath)@mat[referenceRow, ]
  silenceDoGC()
  return(tempMatrix)
})
#gctxMatrix = do.call('cbind', dataList)
gctxMatrix = dataList[[1]]
rm(dataList)
silenceDoGC()

# Data Matrix Trimming
gctxMatrix = gctxMatrix[, sampleDF$inst_id]
silenceDoGC()

# ----- Differential Expression Analysis -----
# Calculate Differential Gene Expression By Condition
cat('Calculating Differential Expression.\n')

# Initiate Multi-Thread
.startTime = date()
currentSPID = SysGetSPID()
referenceColumn = matchedDF$case_ID
caseCount = length(matchedDF$case_ID)
geneCount = length(referenceRow)
registerDoParallel(cores = 120)

# Cleaning Cache
allPath = c(
  .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_fc', idx = 1:caseCount),
  .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_ts', idx = 1:caseCount),
  .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_pv', idx = 1:caseCount),
  .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_fdr', idx = 1:caseCount)
)
silenceRemove(paths = allPath, check = TRUE)

# Parallel Processing
cat('Begin DE Job: ')
tempOutput = foreach(i = 1:caseCount, .inorder = TRUE) %dopar% {
  # Status Update
  if (i %% 500 == 0) {
    cat(i, ' . ', sep = '')
  }
  
  # Obtain Matching-Sample Mapping
  currentTask = matchedDF[i, ]
  trtSamples = sampleDF[.(currentTask$rna_centre, currentTask$trt_type, currentTask$trt_cval), inst_id]
  ctlSamples = sampleDF[.(currentTask$rna_centre, currentTask$ctl_type, currentTask$ctl_cval), inst_id]
  
  # Sample-Size Control
  if (currentTask$trt_orig_sample_size != currentTask$trt_final_sample_size) {
    tempLogical = sample(x = 1:currentTask$trt_orig_sample_size, size = currentTask$trt_orig_sample_size, replace = FALSE)
    tempLogical = (tempLogical %in% 1:currentTask$trt_final_sample_size)
    trtSamples = trtSamples[tempLogical]
    rm(tempLogical)
  }
  
  if (currentTask$ctl_orig_sample_size != currentTask$ctl_final_sample_size) {
    tempLogical = sample(x = 1:currentTask$ctl_orig_sample_size, size = currentTask$ctl_orig_sample_size, replace = FALSE)
    tempLogical = (tempLogical %in% 1:currentTask$ctl_final_sample_size)
    ctlSamples = ctlSamples[tempLogical]
    rm(tempLogical)
  }
  
  # Obtain Data
  trtData = gctxMatrix[, trtSamples]
  ctlData = gctxMatrix[, ctlSamples]
  
  # Calculate: Fold Change
  tempFC = rowMeans(trtData) - rowMeans(ctlData)
  
  # Calculate: T-statistic, P-value, False Discovery Rate
  tempResult = lapply(1:geneCount, function(j) {
    tryCatch(
      {
        tempTest = t.test(x = trtData[j, ], y = ctlData[j, ], alternative = 'two.sided', paired = FALSE, var.equal = TRUE)
        tempVector = c(
          TS = as.numeric(tempTest$statistic),
          PV = tempTest$p.value
        )
        return(tempVector)
      },
      error = function(e) {
        return(c(TS = NA, PV = NA))
      }
    )
  })
  tempTS = unlist(lapply(tempResult, function(j) {return(j['TS'])}))
  tempPV = unlist(lapply(tempResult, function(j) {return(j['PV'])}))
  tempFDR = p.adjust(tempPV, method = 'fdr')
  
  # Save File
  tempPath = .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_fc', idx = i)
  XZSaveRDS(obj = tempFC, file = tempPath)
  
  tempPath = .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_ts', idx = i)
  XZSaveRDS(obj = tempTS, file = tempPath)
  
  tempPath = .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_pv', idx = i)
  XZSaveRDS(obj = tempPV, file = tempPath)
  
  tempPath = .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_fdr', idx = i)
  XZSaveRDS(obj = tempFDR, file = tempPath)
  
  # Variable Cleanup
  rm(currentTask, trtSamples, ctlSamples)
  rm(trtData, ctlData)
  rm(tempResult, tempFC, tempTS, tempPV, tempFDR)
  silenceDoGC()
  return(NULL)
}

# stop the multi-threads
registerDoSEQ()

#clean-up the variable 
rm(matchedDF, sampleDF, gctxMatrix, allPath, tempOutput)
silenceDoGC()
cat('\clean-up is ok.\n')

#Fold Change
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

tempMatrix = foreach(i = 1:caseCount, .inorder = TRUE, .combine = cbind, .maxcombine = 1000) %do% {
  tempPath = .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_ts', idx = i)
  return(readRDS(tempPath))
}
colnames(tempMatrix) = referenceColumn
rownames(tempMatrix) = referenceRow

tempPath = sprintf('%sts.matrix.RDS.XZ', outDir)
XZSaveRDS(obj = tempMatrix, file = tempPath)
rm(tempMatrix)
silenceDoGC()
silenceRemove(paths = .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_ts', idx = 1:caseCount))

tempMatrix = foreach(i = 1:caseCount, .inorder = TRUE, .combine = cbind, .maxcombine = 1000) %do% {
  tempPath = .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_pv', idx = i)
  return(readRDS(tempPath))
}
colnames(tempMatrix) = referenceColumn
rownames(tempMatrix) = referenceRow

tempPath = sprintf('%spv.matrix.RDS.XZ', outDir)
XZSaveRDS(obj = tempMatrix, file = tempPath)
rm(tempMatrix)
silenceDoGC()
silenceRemove(paths = .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_pv', idx = 1:caseCount))

tempMatrix = foreach(i = 1:caseCount, .inorder = TRUE, .combine = cbind, .maxcombine = 1000) %do% {
  tempPath = .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_fdr', idx = i)
  return(readRDS(tempPath))
}
colnames(tempMatrix) = referenceColumn
rownames(tempMatrix) = referenceRow

tempPath = sprintf('%sfdr.matrix.RDS.XZ', outDir)
XZSaveRDS(obj = tempMatrix, file = tempPath)
rm(tempMatrix)
silenceDoGC()
silenceRemove(paths = .Path001(eachThreadDir = ramDir, spid = currentSPID, type = 'lincs_fdr', idx = 1:caseCount))

# Directory Cleanup
silenceRemove(ramDir)

# Print Timestamp
cat(sprintf('START TIME: %s\n', .startTime))
cat(sprintf('END TIME: %s\n\n', date()))


#START TIME: Tue Apr 19 17:46:40 2022
#END TIME: Thu Apr 21 17:31:02 2022
XZSaveRDS(obj = finalSampleMetadata, file = tempPath)

# Matched Metadata
tempPath = sprintf('%scondition.metadata.RDS.XZ', outDir)
XZSaveRDS(obj = finalMatchedMetadata, file = tempPath)

# cat(sprintf('START TIME: %s\n', .startTime))
# START TIME: Fri Apr  1 08:08:28 2022
# cat(sprintf('END TIME: %s\n\n', date()))
# END TIME: Fri Apr  1 23:23:49 2022
################################################################################
################################################################################


###stopCluster(cl)
##
#fc.matrix <- readRDS("/home/xiay/1_home_xiay/Project/PD1/R/35_LINCS_2020/Data/LINCS/Level3/Compound_Data/fc.matrix.RDS.XZ")
#fc.compound_pert.gtc<- new("GCT", mat=fc.matrix)
#write_gctx(fc.compound_pert.gtc, "../../35_LINCS_2020/Data/LINCS/Level3/Compound_Data/fc.compound_pert.gtcx")

setwd("/home/xiay/1_home_xiay/Project/PD1/R/37_LINCS_2020_Results_New/")
metadata.condition <- readRDS("/home/xiay/1_home_xiay/Project/PD1/R/35_LINCS_2020/Data/LINCS/Level3/Compound_Data/condition.metadata.RDS.XZ")


data_trt_cp <- list()
data_trt_cp.df <- list()
fgsea.sam.trt_cp <- list()
fgsea.res.trt_cp <- list()
trt_cp_number.group <- list()

trt_cp_number <- which(metadata.condition$trt_type == "trt_cp")

for(i in 1:54){trt_cp_number.group[[i]] <- trt_cp_number[((i-1)*10000+1):(i*10000)]}
trt_cp_number.group[[55]] <- trt_cp_number[540001:551437]

example_2020 <- parse_gctx("/home/xiay/1_home_xiay/Project/PD1/R/35_LINCS_2020/Data/LINCS/Level3/Compound_Data/fc.compound_pert.gtcx_n551437x12328.gctx",
                           rid=1:12328, cid=1:10)
example_2020 <- example_2020@mat %>% as.data.frame() 
example_2020 <- example_2020 %>% dplyr::mutate(gene_id=rownames(example_2020 ))
geneinfo_beta <- read.delim("/home/xiay/1_home_xiay/Project/PD1/LINCS/3_LINCS_2020/level3/geneinfo_beta.txt")
example_2020$gene_id <- as.integer(example_2020$gene_id)
example_2020 <- dplyr::left_join(example_2020,geneinfo_beta,by= "gene_id")
#there are two genes which the symbol are identical. This will cause the problem, so modified one of them into XXX.2
example_2020[1195,]$gene_symbol <- "MIA2.2" 

saveRDS(example_2020,"example_2020_new.RDS")
save.image("LINCS.2020.cp.Standard.RData")

#################################
#Use fgsea tu perform GSEA 
options(
  stringsAsFactors = FALSE,
  warn = 0
)
library(fgsea)
library(cmapR)
library(tidyverse)

setwd("/home/xiay/1_home_xiay/Project/PD1/R/37_LINCS_2020_Results_New/cp/")
load("./LINCS.2020.cp.Standard.RData")

#library(tictoc)
library(furrr)
library(future)
plan(multisession, workers = 120)

.startTime = date()

for(j in c(1:55)){
  
  load("./LINCS.2020.cp.Standard.RData")
  
  fun_cmap_fgsea <- function(x){
    x <- x %>% unlist() 
    names(x) <-example_2020$gene_symbol
    fgsea.res <- fgsea(pathways = genesets, stats = x,eps= 0.0, minSize  = 5, maxSize  = 500)
    return(fgsea.res)
  }
  genesets <- readRDS("../GSEA_Genesets/super.third.human.pd1.all.RDS")
  data_trt_cp[[j]] <- parse_gctx("/home/xiay/1_home_xiay/Project/PD1/R/35_LINCS_2020/Data/LINCS/Level3/Compound_Data/fc.compound_pert.gtcx_n551437x12328.gctx", 
                                 rid=1:12328, cid=trt_cp_number.group[[j]])
  data_trt_cp.df[[j]] <- as.data.frame(data_trt_cp[[j]]@mat)
  fgsea.res.trt_cp[[j]] <- furrr::future_map(data_trt_cp.df[[j]], ~ fun_cmap_fgsea(.x))
  saveRDS(fgsea.res.trt_cp[[j]],paste("c",j,"_fgsea.res.trt_cp_super.RDS",sep = ""));rm(list=ls());
}

cat(sprintf('START TIME: %s\n', .startTime))
cat(sprintf('END TIME: %s\n\n', date()))


################
#tidy the results

setwd("/home/xiay/1_home_xiay/Project/PD1/R/37_LINCS_2020_Results_New/cp/")
load("/home/xiay/1_home_xiay/Project/PD1/R/37_LINCS_2020_Results_New/cp/LINCS.2020.cp.Standard.RData")
fgsea.res.trt_cp <- list()
for(i in 1:5){fgsea.res.trt_cp[[i]] <- readRDS(paste("/home/xiay/1_home_xiay/Project/PD1/R/37_LINCS_2020_Results_New/cp/c",
                                                     i,"_fgsea.res.trt_cp_super.RDS",sep = ""))}

for(i in 1:8){print(i);
  fun_tidy= function(GSEA_res_list){
    x <- furrr::future_map_dfr(GSEA_res_list, ~ .x[i,])
    return(x)
  }
  fgsea.res.tidy[[i]] <- furrr::future_map(fgsea.res.trt_cp,~ fun_tidy(.x)) %>% Reduce(rbind,.)
}
#furrr::future_map_dfr(fgsea.res.trt_cp[[1]],~ .x[i,]) %>% View()
cl <- makeCluster(8)
library(doParallel)
library(foreach)
library(furrr)
registerDoParallel(cl)
fgsea.res.tidy <- list()
fgsea.res.tidy.list<- list()
plan(multisession, workers =8)


#######
#
.startTime = date()
fgsea.res.tidy.try <- list()
library(doParallel)
registerDoParallel(11)
fgsea.res.tidy.try <- foreach (m = 1:11) %dopar% {
  fgsea.res.trt_cp <- list()
  k=1;for(n in c(((m-1)*5+1):(m*5))){fgsea.res.trt_cp[[k]] <- readRDS(paste(
    "/home/xiay/1_home_xiay/Project/PD1/R/37_LINCS_2020_Results_New/cp/c",n,"_fgsea.res.trt_cp_super.RDS",sep = ""));k=k+1}
  library(doParallel)
  registerDoParallel(8)
  fgsea.res.tidy.list <- foreach(i = 1:8 ) %dopar% {
    fun_tidy= function(GSEA_res_list){
      library(purrr)
      #plan(multisession, workers =8)
      
      x <- purrr::map_df(GSEA_res_list, ~ .x[i,])
      #x <- furrr::future_map_dfr(GSEA_res_list, ~ .x[i,])
      return(x)
    }
    library(tidyverse)
    purrr::map(fgsea.res.trt_cp,~ fun_tidy(.x )) %>% Reduce(rbind,.)
    #fgsea.res.tidy[[i]] %>% Reduce(rbind,fgsea.res.tidy.list)
  } 
  #fgsea.res.tidy.try[[m]] <- fgsea.res.tidy.list
}

cat(sprintf('START TIME: %s\n', .startTime))
cat(sprintf('END TIME: %s\n\n', date()))


list_temp <- list()
for(i in 1:8){list_temp[[i]] <- list()}

for(j in 1:8){
  for(i in 1:11){
    list_temp[[j]][[i]] <-  fgsea.res.tidy.try[[i]][[j]]
  } 
}      

fgsea.clue.order <- list()
for(i in 1:8){
  fgsea.clue.order[[i]] <- Reduce(rbind,list_temp[[i]])
}
fgsea.clue.order_with_id <- list()

for(i in 1:8){fgsea.clue.order[[i]] <- fgsea.clue.order[[i]] %>% dplyr::mutate(id=rownames(fgsea.clue.order[[i]]))}
#saveRDS(fgsea.clue.order_with_id,"fgsea.clue.order_with_id.RDS")

fgsea.clue.screen1 <- purrr::map(fgsea.clue.order,~ dplyr::filter(.x,padj<0.05,NES>0))
fgsea.clue.screen2 <- purrr::map(fgsea.clue.screen1,~  .x$id)
fgsea.clue.screen3 <- Reduce(intersect,fgsea.clue.screen2)



##screen data
setwd("/home/xiay/1_home_xiay/Project/PD1/R/37_LINCS_2020_Results_New/1_cp/")
load("./LINCS.2020.cp.Standard.RData")

list.good <- list()
for(i in 1:8){list.good[[i]] <- fgsea.clue.order[[i]][as.numeric(good.id),]}
good.df <- as.data.frame(matrix(data = NA,nrow =dim(list.good[[1]])[1],ncol = 8));{for(i in 1:8){good.df[,i] <- list.good[[i]][,6]}}
colnames(good.df) <- names(fgsea.clue.order)

good.score <- vector();
####CM-Score
for(i in 1:dim(list.good[[1]])[1]){good.score[i] <- 0.4986+0.0969*good.df[,1][i]+0.0892*good.df[,2][i]
+0.0307*good.df[,3][i]+0.0117*good.df[,4][i]+0.0124*good.df[,6][i]+ 0.0213*good.df[,7][i]}

metadata.condition <- readRDS("/home/xiay/1_home_xiay/Project/PD1/R/35_LINCS_2020/Data/LINCS/Level3/Compound_Data/condition.metadata.RDS.XZ")
good.df.withscore <- data.frame(good.score,good.df)
good.df.withscore <- data.frame(rownames(good.df.withscore),good.df.withscore)
colnames(good.df.withscore)[1] <- "good_id_index"
good.df.withscore.order <- good.df.withscore %>% arrange(desc(good.score))

good.id.order <- good.id[as.numeric(good.df.withscore.order$good_id_index)]
good.id.order.trt_number <- trt_cp_number[as.numeric(good.id.order)]

#metadata.condition <- readRDS("../31_LINCS_Phase1_level3toFC/Data/LINCS/Level3/condition.metadata.RDS.XZ")
sig_info.choose <- metadata.condition[good.id.order.trt_number,]
sig_info.choose <- data.frame(1:dim(list.good[[1]])[1],sig_info.choose)
colnames(sig_info.choose)[1] <- "id"
sig_info.choose <- sig_info.choose %>% mutate(Pert_Score=good.df.withscore.order$good.score)
sig_info.choose.fgsea <- cbind(sig_info.choose,good.df.withscore.order)

#########
cellinfo_beta <- read.delim("./1_home_xiay/Project/PD1/LINCS/3_LINCS_2020/level3/cellinfo_beta.txt")

cell_line <- table(sig_info.choose.fgsea$cell_id) %>% as.data.frame() %>% arrange(-Freq)
colnames(cell_line)[2] <- "Freq"
colnames(cell_line)[1] <- "cell_iname"
cellinfo_beta.tumor <- cellinfo_beta %>% dplyr::filter(cell_type=="tumor")
cellinfo_beta.tumor.rm_unknown <- cellinfo_beta.tumor %>% dplyr::filter(primary_disease!="unknown")

cell_line.tumor.rm_unknown.index <-  which(cellinfo_beta.tumor.rm_unknown$cell_iname %in% cell_line$cell_iname)
cell_line.tumor.rm_unknown <- cellinfo_beta.tumor.rm_unknown[cell_line.tumor.rm_unknown.index,]

table(cellinfo_beta.tumor.rm_unknown$primary_disease) %>% as.data.frame() %>%  dplyr::arrange(-Freq)




primary_disease  <- cellinfo_beta.tumor.rm_unknown$primary_disease %>%  unique()

cell_iname.use <- cellinfo_beta.tumor.rm_unknown$cell_iname


cell_line.use <- cell_line[which(cell_line$cell_iname %in% cell_iname.use),]

#group by cancer type
cell_line.use.list <- purrr::map(primary_disease,~dplyr::filter(cellinfo_beta.tumor.rm_unknown,primary_disease==.) 
                                 %>% .$cell_iname)

sig_info.choose.fgsea.list <- purrr::map(cell_line.use.list,~sig_info.choose.fgsea[which(sig_info.choose.fgsea$cell_id %in% .),])

drug.choose.list <- purrr::map(sig_info.choose.fgsea.list,~.$trt_iname)




##########################################################################################
##########################################################################################
#Screen the small molecular compounds in A549 cell line (or A375 cell line)
sig_info.choose.A549<- sig_info.choose.fgsea %>% dplyr::filter(cell_id == "A549")
drug.choose.detail <- table(sig_info.choose$trt_iname) %>% as.data.frame()%>% arrange(desc(Freq))
drug.choose.detail <- drug.choose.detail %>% mutate(trt_times=NA,trt_cell_kind=NA,
                                                    choose_cell_kind=NA,A549=NA,A549_times=NA,
                                                    pert_score_max=NA,pert_score_detail=NA)

#plan(multisession, workers = 20)

colnames(drug.choose.detail)[1] <- "compound"
drug.choose.detail$trt_times <- purrr::map(drug.choose.detail$compound,~ sum(metadata.condition$trt_iname==.x))  %>% unlist()
drug.choose.detail$trt_cell_kind <- purrr::map(
  drug.choose.detail$compound,~ dplyr::filter(metadata.condition,trt_iname==.x) 
  %>% select(cell_id) %>% table() %>% as.data.frame() %>% dim() %>% .[1])


drug.choose.detail$choose_cell_kind <- purrr::map(
  drug.choose.detail$compound,~ dplyr::filter(sig_info.choose,trt_iname==.x) 
  %>% select(cell_id) %>% table() %>% as.data.frame() %>% dim() %>% .[1])
drug.choose.detail$choose_cell_kind <- drug.choose.detail$choose_cell_kind %>% unlist()



#plan(multisession, workers = 32)

for(i in c(6)){
  temp <- furrr::future_map(
    drug.choose.detail$compound,~ dplyr::filter(sig_info.choose,trt_iname==.x,cell_id==colnames(drug.choose.detail)[i])
    %>% dim() %>% .[1])
  drug.choose.detail[,i] <- unlist(temp)
}

A549_times <- furrr::future_map(
  drug.choose.detail$compound,~ dplyr::filter(metadata.condition,trt_iname==.x,cell_id=="A549")
  %>% dim() %>% .[1])
drug.choose.detail$A549_times<- unlist(A549_times)
drug.choose.detail$A549_times <-drug.choose.detail$A549_times %>%  unlist()





drug.choose.detail$pert_score_max <-purrr::map(drug.choose.detail$compound,~ dplyr::filter(sig_info.choose,trt_iname==.x) 
                                               %>% select(Pert_Score)%>% max() %>% round(2) %>% unlist() )
drug.choose.detail$pert_score_max<- drug.choose.detail$pert_score_max %>% unlist()

drug.choose.detail$pert_score_detail <-purrr::map(drug.choose.detail$compound,~ dplyr::filter(sig_info.choose,trt_iname==.x) %>% 
                                                    select(Pert_Score)%>% round(2) %>% unlist(.) %>% as.data.frame() %>% .[,1]  )

fun_paste <- function(x,y){
  paste(x,y,sep="_")
}

drug.choose.detail$pert_score_detail <- purrr::map(drug.choose.detail$pert_score_detail,~ Reduce(fun_paste,.))

drug.choose.detail <- purrr::map_df(drug.choose.detail,~ unlist(.)) 
drug.choose.detail[is.na( drug.choose.detail)] <-  -666


#.libPaths(c("/usr/local/lib/R/site-library","/usr/local/lib/R/library","./Rpackage/","/home/xiay/1_home_xiay/Software/R/site-library/"))
library(openxlsx)
#save in excel file-type
openxlsx::write.xlsx(x = drug.choose.detail , file = "./results_with_CM_Drug.xlsx",
                     sheetName = "screenResult", rownames = FALSE)
