rm(list=ls())

library(DESeq2)
library(biomaRt)
library(combinat)
library(dplyr)

maximal_humans_to_consider_per_condition <- 6 


out_base = "/Users/tstoeger/Dropbox/aging_map_paper/dynamic/tstoeger/200609_gtex_f"
pfus_of_interest = c(0)



sort_each_half <- function(v){
  
  v <- as.numeric(unlist(v))
  samples = length(v)
  half = ceiling(samples/2)
  
  first <- sort(v[1:half])
  last <- sort(v[(half+1):samples])
  result <- c(first, last)
  return(result)
}


get_combinations <- function(replicates){
  permutations <- permn(replicates) 
  sorted_combos <- lapply(permutations, sort_each_half)
  sorted_combos <- unique(sorted_combos)
  
  return(sorted_combos)
}

get_half <- function(combo, first_or_last){
  combo <- unlist(combo)
  samples <- length(combo)
  half <- ceiling(samples/2) 
  
  if (first_or_last=='first'){
    combo <- combo[1:half]
  } else if (first_or_last=='last'){
    combo <- combo[(half+1):samples]
  } else{
    stop('first_or_last must either be first or last')
  }
  return(combo)
}

# Prepare output folde in case it was missing
if (dir.exists(out_base)==FALSE){
  dir.create(out_base, recursive = TRUE)            
}

#load csv of overview metrics
overview <- read.csv(file.path(out_base, "sample_meta.csv"))
groups <- overview[,c('orig_sample_name', 'tissue', 'pfu', 'age', 'human_id')]
groups <- groups[which(groups$pfu %in% pfus_of_interest), ]

#select counts files for analysis
files <- paste0(groups$orig_sample_name,".counts")

#go to directory with htseq counts files and normalize subsets corresponding to cell types in celllist
setwd(file.path(out_base, 'counts/'))

bm <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
AnnotationBM <- getBM(mart=bm, attributes=c('ensembl_gene_id','gene_biotype','external_gene_name','description'))

# Initialize output directory
out_dir = file.path(out_base, "DE/Flu/")
if (dir.exists(out_dir)==FALSE){
  dir.create(out_dir, recursive = TRUE)            
}


max_combo = 0




## cycle through comparisons
for (tissue in unique(groups$tissue)){
  
  stable <- groups[groups$tissue == tissue,]
  
  for (pfu in unique(stable$pfu)){
    
    ttable <- stable[stable$pfu == pfu,]
    ttable <- ttable[,c('orig_sample_name','age', 'human_id')]
    
    for (age in unique(ttable$age)){
      utable <- ttable[ttable$age == age,]
      
      
      
      
      number_of_replicates = length(utable$human_id)
      stopifnot(number_of_replicates==length(unique(utable$human_id)))
      
      if (number_of_replicates>=4){
        
        if (dim(utable)[1]>maximal_humans_to_consider_per_condition){
          utable <- sample_n(utable, maximal_humans_to_consider_per_condition,replace = FALSE)
        }
        
        
        
        combos <- get_combinations(unique(utable$human_id))
        
        
        for (combo in combos){
          
          vtable <- utable
          
          firsts <- get_half(combo, 'first')
          lasts <- get_half(combo, 'last')
          
          vtable$condition[vtable$human_id %in% firsts] <- 1
          vtable$condition[vtable$human_id %in% lasts] <- 2
          
          if (any(is.na(vtable$condition))){
            stop('At least one mouse could not be matched.')
          }
          
          vtable <- vtable[c('orig_sample_name', 'condition')]
          colnames(vtable) <- c("samples","condition")
          rownames(vtable) <- vtable$samples
          
          condition <- vtable$condition   #this is the groups file for deseq
          
          print(
            paste0(
              "comparisons for ",tissue," flu level ", pfu," in age ", age, " in mice ",
              paste(firsts,collapse = " "), ' against ', paste(lasts,collapse = " "))
          )
          
          
          tmp <- files[match(paste0(rownames(vtable),".counts"),files)]
          
          datalist <- lapply(tmp, function(x){read.delim(x, header = F,row.names = 1)})
          
          rawdata <- Reduce(function(x,y) {data.frame(x,y)}, datalist)
          colnames(rawdata) <- rownames(vtable) #this is the final input counts file
          
          safety_check_colnames <- colnames((rawdata))   # the first row appears to carry sample names -> remove safely
          safety_check_first_row <- names(unlist(rawdata[1,]))
          stopifnot(all(safety_check_colnames==safety_check_first_row))
          
          
          rawdata <- rawdata[c(2:(nrow(rawdata))),]
          
          # indx <- sapply(rawdata, is.factor)    # convert to numbers
          # rawdata[indx] <- lapply(rawdata[indx], function(x) as.numeric(as.character(x)))
          # 
          # 
          # rawdata <- lapply(rawdata, function(x) as.numeric(as.character(x)))
          # 
          
          # b <- rawdata
          # 
          # for (x in colnames(rawdata)){
          #   print(x)
          #   
          # }
          # 
          indx = colnames(rawdata)
          rawdata[indx] <- lapply(rawdata[indx], function(x) as.numeric(as.character(x)))
          
          
          
          # a <- DESeqDataSetFromMatrix(countData=b, colData=vtable, design=~condition)
          
          
          
          
          dds <- DESeqDataSetFromMatrix(countData=rawdata, colData=vtable, design=~condition)
          dds <- dds[ rowSums(counts(dds)) > 0, ]    # changed in the any detected versions from 1 to 0
          rld <- rlog(dds, blind=TRUE)
          
          # DESEq analysis
          dds<-DESeq(dds)
          res<-results(dds,alpha=0.05) #change cutoffs?
          summary(res)
          ensembl_ids <- rownames(assay(rld))
          stopifnot(rownames(res)==rownames(rawdata[ensembl_ids,]))
          res1 <- cbind(as.data.frame(row.names(rawdata[ensembl_ids,])),assay(rld),res)
          
          # Determine differentially experssed genes
          up <- res1$padj < 0.05 & res1$log2FoldChange > 0
          flag <- rep(0, nrow(res1))
          flag[up] <- 1
          res1$up <- flag
          dn <- res1$padj < 0.05 & res1$log2FoldChange < 0
          flag1 <- rep(0, nrow(res1))
          flag1[dn] <- 1
          res1$dn <- flag1
          res1 <- res1[order(-res1$up, -res1$dn, res1$padj),]
          colnames(res1)[1] <- "Symbol"
          
          # Export
          file_name = paste0(tissue,"_pfu_",pfu,"_age_",age, '_first_', paste(firsts,collapse = "-"), "_DE.csv")
          p = file.path(out_dir, file_name)
          write.csv(res1,p,row.names = T)
          
        }
        
      }
    }
  }
}

savehistory(file.path(out_base, "DE/datalog.txt"))
