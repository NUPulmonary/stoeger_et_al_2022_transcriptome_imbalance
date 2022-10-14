rm(list=ls())

library(DESeq2)
library(biomaRt)
library(combinat)

out_base = "/Users/tstoeger/Dropbox/aging_map_paper/dynamic/tstoeger/200129_inner_bootstrap"
tissues_of_interest = c(
  'Lung', 'Blood', 'Brain', 'Cerebellum', 'Heart', 'Kidney', 'GutEP', 
  'Liver', 'Adrenal', 'Esophagus', 'Stomach', 'SI', 'BAT', 'WAT', 'LI', 'Skin', 'MuscSat')
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
  # 
  # combos = list()
  # 
  # samples <- length(replicates)
  # half <- ceiling(samples/2)
  # 
  # for (j in 1:length(sorted_combos)){
  #   v = unlist(sorted_combos[j])
  #   first <- v[1:half]
  #   for(jj in 1:length(sorted_combos)){
  #     v = unlist(sorted_combos[jj])
  #     last <- v[(half+1):samples]
  #     
  #     
  #     could_export <- FALSE
  #     if (samples%%2==0){
  #       if ((all(first==last)) && (j<jj)){
  #         could_export <- TRUE
  #       }
  #     }else if (samples%%2==1){
  #       if (j==jj){
  #         could_export <- TRUE
  #       }
  #     }else{
  #       stop('Sample number is not anticipated by code.')
  #     }
  #     
  #     
  #     if (could_export){
  #       to_insert <- sorted_combos[j]
  #       combos <- append(combos, to_insert)
  #     }
  #   }
  # }
  # 
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
groups <- overview[,c('orig_sample_name', 'tissue', 'pfu', 'age', 'mouse_id', 'recommend_to_discard_181022')] 
groups <- groups[-c(which(groups$recommend_to_discard_181022==1 | groups$age==.5)),]
groups <- groups[,c('orig_sample_name', 'tissue', 'pfu', 'age', 'mouse_id')]
groups <- groups[which(groups$tissue %in% tissues_of_interest), ]
groups <- groups[which(groups$pfu %in% pfus_of_interest), ]

#select counts files for analysis
files <- paste0(groups$orig_sample_name,".htseq.counts")

#go to directory with htseq counts files and normalize subsets corresponding to cell types in celllist
setwd(file.path(out_base, 'mockquestmon2_counts/'))

bm <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
AnnotationBM <- getBM(mart=bm, attributes=c('ensembl_gene_id','gene_biotype','external_gene_name','description'))

# Initialize output directory
out_dir = file.path(out_base, "DE/Flu/")
if (dir.exists(out_dir)==FALSE){
  dir.create(out_dir, recursive = TRUE)            
}


## cycle through comparisons
for (tissue in unique(groups$tissue)){
  
  stable <- groups[groups$tissue == tissue,]

  for (pfu in unique(stable$pfu)){
    
    ttable <- stable[stable$pfu == pfu,]
    ttable <- ttable[,c('orig_sample_name','age', 'mouse_id')]
    
    for (age in unique(ttable$age)){
      utable <- ttable[ttable$age == age,]

      
      
      
      number_of_replicates = length(utable$mouse_id)
      stopifnot(number_of_replicates==length(unique(utable$mouse_id)))
          
      if (number_of_replicates>=2){
        combos <- get_combinations(unique(utable$mouse_id))
        
        for (combo in combos){
    
          vtable <- utable
          
          firsts <- get_half(combo, 'first')
          lasts <- get_half(combo, 'last')
          
          vtable$condition[vtable$mouse_id %in% firsts] <- 1
          vtable$condition[vtable$mouse_id %in% lasts] <- 2
          
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
          
          
          tmp <- files[match(paste0(rownames(vtable),".htseq.counts"),files)]
          datalist <- lapply(tmp, function(x){read.delim(x, header = F,row.names = 1)})    
          rawdata <- Reduce(function(x,y) {data.frame(x,y)}, datalist)
          colnames(rawdata) <- rownames(vtable) #this is the final input counts file
          
          safety_check_colnames <- colnames((rawdata))   # the first row appears to carry sample names -> remove safely
          safety_check_first_row <- names(unlist(rawdata[1,])) 
          stopifnot(all(safety_check_colnames==safety_check_first_row))
          rawdata <- rawdata[c(2:(nrow(rawdata))),]
          
          indx <- sapply(rawdata, is.factor)    # convert to numbers
          rawdata[indx] <- lapply(rawdata[indx], function(x) as.numeric(as.character(x)))
          
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
