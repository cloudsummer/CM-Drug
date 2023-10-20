

################################################################################
A375 Cell Line (Using Level5 data)
################################################################################
rm(list = ls())
setwd("workingDirectory/")
CM_Drug.condition <- fread(input = "./pertInfo.Data2.A375", sep = '\t', header = TRUE, data.table = FALSE)

data_trt_cp <- list()
data_trt_cp.df <- list()
fgsea.sam.trt_cp <- list()
fgsea.res.trt_cp <- list()
trt_cp_number.group <- list()

trt_cp_number <- which(CM_Drug.condition$pert_type == "trt_cp")

trt_cp_number.group[[1]] <- trt_cp_number[1:length(trt_cp_number)]

template <- parse_gctx("./Data_A375.2.gctx",
                       rid=1:12328, cid=1:10)@mat %>% as.data.frame()

template <- template %>% dplyr::mutate(gene_id=rownames(template ))

gene_df <- read.delim("./gene_df.2")

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
setwd("workingDirectory/")
resultDir="./result/cp/"
dir.create(resultDir, recursive = TRUE, showWarnings = FALSE)
load("./template.RData")

#library(tictoc)
library(furrr)
library(future)
plan(multisession, workers = 120)

.startTime = date()

for(j in c(1)){
  setwd("workingDirectory/")
  load("./template.RData")
  
  fun_cmap_fgsea <- function(x){
    x <- x %>% unlist() 
    names(x) <-template$gene_symbol
    fgsea.res <- fgsea(pathways = genesets, stats = x,eps= 0.0, minSize  = 5, maxSize  = 500)
    return(fgsea.res)
  }
  genesets <- readRDS("super.third.human.pd1.all.RDS")
  data_trt_cp[[j]] <- parse_gctx("./Data_A375.2.gctx", 
                                 rid=1:12328, cid=trt_cp_number.group[[j]])
  data_trt_cp.df[[j]] <- as.data.frame(data_trt_cp[[j]]@mat)
  fgsea.res.trt_cp[[j]] <- furrr::future_map(data_trt_cp.df[[j]], ~ fun_cmap_fgsea(.x))
  setwd("./result/cp/");saveRDS(fgsea.res.trt_cp[[j]],paste("c",j,"_fgsea.res.trt_cp_super.RDS",sep = ""));rm(list=ls());
}

cat(sprintf('START TIME: %s\n', .startTime))
cat(sprintf('END TIME: %s\n\n', date()))

################################################################################
#tidy the results

setwd("workingDirectory/result/cp/")
load("workingDirectory/template.RData")

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

setwd("workingDirectory/")

fgsea.clue.order <- readRDS("./result/cp/fgsea.clue.order_with_id.RDS")

#screening strategy1
fgsea.clue.screen1 <- purrr::map(fgsea.clue.order,~ dplyr::filter(.x,padj<0.05,NES>0))
fgsea.clue.screen2 <- purrr::map(fgsea.clue.screen1,~  .x$id)
fgsea.clue.screen3 <- Reduce(intersect,fgsea.clue.screen2)
good.id <- fgsea.clue.screen3

#screening strategy2
#Considering that certain drugs significantly enhance specific pathways but show slightly insufficient significance in one or two other pathways, to prevent the oversight of such drugs, we have implemented an alternative threshold screening allowing for some flexibility. Specifically, we allow an adjusted p-value between 0.05 and 0.2 for enrichment results in two pathways, while the adjusted p-value for other enriched pathways must be less than 0.05.
# fgsea.clue.screen1 <- purrr::map(fgsea.clue.order,~ dplyr::filter(.x,padj<0.2,NES>0))
# fgsea.clue.screen2 <- purrr::map(fgsea.clue.screen1,~  .x$id)
# fgsea.clue.screen3 <- Reduce(intersect,fgsea.clue.screen2)
# fgsea.clue.screen4 <- purrr::map(fgsea.clue.order,~ .x[as.integer(fgsea.clue.screen3),])
# fgsea.clue.screen5 <- purrr::map(fgsea.clue.screen4 ,~ .x[,3])
# fgsea.clue.screen6 <- Reduce(cbind,fgsea.clue.screen5)
# fgsea.clue.screen6 <- cbind(fgsea.clue.screen6,fgsea.clue.screen4[[1]]$id)
# fgsea.clue.screen6 <- fgsea.clue.screen6 %>% as.data.frame()
# fgsea.clue.screen6.row <- purrr::map(as.data.frame(t(fgsea.clue.screen6[,1:8])), ~ .x)
# fgsea.clue.screen6.row <- purrr::map(fgsea.clue.screen6.row,~ as.numeric(.))
# fgsea.clue.screen6.count <- purrr::map(fgsea.clue.screen6.row,function(x){count <- 0;
# for(i in 1:8){if(x[i]<0.05){count <- count+1}};return(count);count <- 0})
# unlist(fgsea.clue.screen6.count) %>% as.data.frame() -> temp
# colnames(fgsea.clue.screen6)[1:8] <- paste(colnames(fgsea.clue.screen6)[1],1:8,sep = "");
# colnames(fgsea.clue.screen6)[9] <- "id"
# fgsea.clue.screen6 <- fgsea.clue.screen6  %>% mutate(num_of_0.05=temp$.)
# fgsea.clue.screen.end <- dplyr::filter(fgsea.clue.screen6,num_of_0.05>=7)
# fgsea.clue.screen.id <- fgsea.clue.screen.end$id
# good.id <- fgsea.clue.screen.id

##screen data
load("./template.RData")

list.good <- list()
for(i in 1:8){list.good[[i]] <- fgsea.clue.order[[i]][as.numeric(good.id),]}
good.df <- as.data.frame(matrix(data = NA,nrow =dim(list.good[[1]])[1],ncol = 8));{for(i in 1:8){good.df[,i] <- list.good[[i]][,6]}}
colnames(good.df) <- names(fgsea.clue.order)

CM_Score <- vector();

################################################################################
#CM-Score

for(i in 1:dim(list.good[[1]])[1]){CM_Score[i] <- 0.4986+0.0969*good.df[,1][i]+0.0892*good.df[,2][i]+0.0307*good.df[,3][i]+0.0117*good.df[,4][i]+0.0124*good.df[,6][i]+ 0.0213*good.df[,7][i]}

CM_Drug.condition <- fread(input = "./pertInfo.Data2.A375", sep = '\t', header = TRUE, data.table = FALSE)

good.df.withscore <- data.frame(CM_Score,good.df)

good.df.withscore <- data.frame(rownames(good.df.withscore),good.df.withscore)

colnames(good.df.withscore)[1] <- "good_id_index"

good.df.withscore.order <- good.df.withscore %>% arrange(desc(CM_Score))

good.id.order <- good.id[as.numeric(good.df.withscore.order$good_id_index)]

good.id.order.trt_number <- trt_cp_number[as.numeric(good.id.order)]

sig_info.choose <- CM_Drug.condition[good.id.order.trt_number,]

sig_info.choose <- data.frame(1:dim(list.good[[1]])[1],sig_info.choose)

colnames(sig_info.choose)[1] <- "id"

finalScreen <- cbind(sig_info.choose,good.df.withscore.order) %>% dplyr::select(3,12:20)

################################################################################
#save in excel file-type
library(openxlsx)
openxlsx::write.xlsx(x = finalScreen , file = "./result/cp/results_with_CM_Drug.2.A375.xlsx",
                     sheetName = "screenResult", rownames = FALSE)
