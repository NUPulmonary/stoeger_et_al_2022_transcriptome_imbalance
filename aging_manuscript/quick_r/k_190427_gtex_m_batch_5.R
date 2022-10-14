# This script mimics the input for DE by pulmonary medicine, as used inhouse
# it still keeps pfu (set to 0)
# order of age and pfu are swapped
# rlog and output dependent on it were removed (after separately testing if same output on pAdj and fold change)
# run on different computer than inital analysis on mouse (iMac rather than macbook, with slightly different versions - after verifying practically 
#   identical output (ratio of fold change and padj is ~1.000change)


rm(list=ls())
library(DESeq2)
library(biomaRt)

bm <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
AnnotationBM <- getBM(mart=bm, attributes=c('ensembl_gene_id','gene_biotype','external_gene_name','description'))

chosen_batch <- 5

out_base = "/Users/tstoeger/Dropbox/aging_map_paper/dynamic/tstoeger/190427_gtex_m/"
if (dir.exists(out_base)==FALSE){
  dir.create(out_base, recursive = TRUE)         
}




#load csv of overview metrics
overview <- read.csv(paste0(out_base, "sample_meta.csv"))


# groups <- overview[,c(2,4,5,7,ncol(overview))]  # now the last column contains a distinct definition
groups <- overview[,c(1,2,3,4,ncol(overview))]  # now the last column contains a distinct definition

groups <- groups[c(which(groups$batch==chosen_batch)),]


# groups <- groups[-c(which(groups$recommend_to_discard==1)),]
# groups <- groups[-c(which(groups$recommend_to_discard_181022==1 | groups$age==.5)),]
#groups$orig_sample_name[which(duplicated(groups$orig_sample_name) == T)]#UMRNA here that we dont need
#groups <- groups[-c(which(duplicated(groups$orig_sample_name) == T)),]



groups <- droplevels(groups[,-ncol(groups)])

#select counts files for analysis (remove anything not in overview csv and with reccomend to discard 1 and age = .5)
# files <- paste0(groups$orig_sample_name,".counts")
files <- paste0(groups$run_name,".counts")


#go to directory with htseq counts files and normalize subsets corresponding to cell types in celllist
setwd(paste0(out_base, 'counts/'))


for (i in unique(groups$tissue)){
  print(i)
}


# #doing for all of same age within a tissue first
# for (i in unique(groups$tissue)){ #i="AM"
#   stable <- groups[groups$tissue == i,]
#   ages <- unique(stable$age)
#   for (j in ages){ #j=9
#     ttable <- stable[stable$age == j,]
#     ttable <- ttable[,-c(2,4)]
#     if (length(unique(ttable$pfu)) >=2){
#       combos <- combn(unique(ttable$pfu),2,simplify = F)
#     } else {combos = NULL}
#     if(length(unique(ttable$pfu)) != 1 & length(combos) >= 1){
#       for (k in combos) {
#         utable <- ttable[which(ttable$pfu %in% unlist(k)),]#k=combos[1]
#         colnames(utable) <- c("samples","condition")
#         rownames(utable) <- utable$samples
#         utable <- utable[-1]
#         utable$condition <- factor(utable$condition)
#         condition <- utable$condition#this is the groups file for deseq
#         
#         print(paste0("comparisons for ",i," aged ", j," across flu levels ", paste(k,collapse = " ")))
#         
#         tmp <- files[match(paste0(rownames(utable),".htseq.counts"),files)]
#         datalist <- lapply(tmp, function(x){read.delim(x, header = F,row.names = 1)})
#         rawdata <- Reduce(function(x,y) {data.frame(x,y)}, datalist)
#         rawdata <- rawdata[c(1:(nrow(rawdata))),] #exclude __no feature,__ambiguous,etc from normalization
#         
#         
#         ### TS: INSERTSTART: CONVERT TO NUMBERS
#         indx <- sapply(rawdata, is.factor)
#         rawdata[indx] <- lapply(rawdata[indx], function(x) as.numeric(as.character(x)))
#         ### TS: INSERT END: CONVERT TO NUMBERS
#         
#         colnames(rawdata) <- rownames(utable) #this is the final input counts file
#         
#         dds <- DESeqDataSetFromMatrix(countData=rawdata, colData=utable, design=~condition)
#         dds <- dds[ rowSums(counts(dds)) > 0, ]           # changed in the any detected versions from 1 to 0
#         rld <- rlog(dds, blind=TRUE)
#         
#         # DESEq analysis
#         dds<-DESeq(dds)
#         res<-results(dds,alpha=0.05) #change cutoffs?
#         summary(res)
#         
#         ensembl_ids <- rownames(assay(rld))
#         stopifnot(rownames(res)==rownames(rawdata[ensembl_ids,]))
#         
#         res1 <- cbind(as.data.frame(row.names(rawdata[ensembl_ids,])),assay(rld),res)
#         #iv <- match(rownames(res1),AnnotationBM$ensembl_gene_id)
#         #res1$gene <- AnnotationBM[iv,'external_gene_name']
#         up <- res1$padj < 0.05 & res1$log2FoldChange > 0 #cutoffs to mess with
#         flag <- rep(0, nrow(res1))
#         flag[up] <- 1
#         res1$up <- flag
#         dn <- res1$padj < 0.05 & res1$log2FoldChange < 0 #also cutoffs to mess with
#         flag1 <- rep(0, nrow(res1))
#         flag1[dn] <- 1
#         res1$dn <- flag1
#         #res1$biotype <- AnnotationBM[iv,'gene_biotype']
#         #res1$description <- AnnotationBM[iv,'description']
#         res1 <- res1[order(-res1$up, -res1$dn, res1$padj),]
#         colnames(res1)[1] <- "Symbol"
#         
#         
#         
#         
#         out_dir = paste0(out_base, "/DE/Age/")
#         if (dir.exists(out_dir)==FALSE){
#           dir.create(out_dir, recursive = TRUE)            
#         }
#         p = paste0(out_dir,i,"_Age_",j,"_pfu_",paste(k,collapse = " "),"_DE.csv")
#         
#         write.csv(res1,p,row.names = T)
#         
#         
#         
#         
#         
#       } 
#     } else {print(paste(i," Age ",j," has only one flu condition"))}
#   }
# }




##now for constant flu condition
for (i in unique(groups$tissue)){ #i="AM"
  stable <- groups[groups$tissue == i,]
  flus <- unique(stable$pfu)
  for (j in flus){ #j=0
    ttable <- stable[stable$pfu == j,]
    # ttable <- ttable[,-c(2,3)]
    ttable <- ttable[,-c(2,4)]
    if (length(unique(ttable$age)) >=2){
      combos <- combn(unique(ttable$age),2,simplify = F)
    } else {combos = NULL}
    if(length(unique(ttable$age)) != 1 & length(combos) >= 1){
      for (k in combos) {
        utable <- ttable[which(ttable$age %in% unlist(k)),]#k=combos[1]
        colnames(utable) <- c("samples","condition")
        rownames(utable) <- utable$samples
        utable <- utable[-1]
        utable$condition <- factor(utable$condition)
        condition <- utable$condition#this is the groups file for deseq

        print(paste0("comparisons for ",i," flu level ", j," across ages ", paste(k,collapse = " ")))

        tmp <- files[match(paste0(rownames(utable),".counts"),files)]

        datalist <- lapply(tmp, function(x){read.delim(x, header = F,row.names = 1)})    # original

        rawdata <- Reduce(function(x,y) {data.frame(x,y)}, datalist)
        ## ORIGINAL # rawdata <- rawdata[c(1:(nrow(rawdata)-5)),] #exclude __no feature,__ambiguous,etc from normalization
        # rawdata <- rawdata[c(2:(nrow(rawdata)-5)),] # TS: also skip first row which would be character of sample labels
        rawdata <- rawdata[c(2:(nrow(rawdata))),] # TS: also skip first row which would be character of sample labels, lo longer skip last ones as they no longer seem in input files

        ### TS: INSERTSTART: CONVERT TO NUMBERS
        indx <- sapply(rawdata, is.factor)
        rawdata[indx] <- lapply(rawdata[indx], function(x) as.numeric(as.character(x)))
        ### TS: INSERT END: CONVERT TO NUMBERS

        colnames(rawdata) <- rownames(utable) #this is the final input counts file

        dds <- DESeqDataSetFromMatrix(countData=rawdata, colData=utable, design=~condition)
        dds <- dds[ rowSums(counts(dds)) > 0, ]            # changed in the any detected versions from 1 to 0
        # TS: 190427 #  rld <- rlog(dds, blind=TRUE)

        # DESEq analysis
        dds<-DESeq(dds)
        res<-results(dds,alpha=0.05) #change cutoffs?
        summary(res)

        #  TS: 190427 # ensembl_ids <- rownames(assay(rld))
        
        ensembl_ids <- rownames(res)
        
        stopifnot(rownames(res)==rownames(rawdata[ensembl_ids,]))

        # TS: 190427 #  res1 <- cbind(as.data.frame(row.names(rawdata[ensembl_ids,])),assay(rld),res)
        res1 <- cbind(as.data.frame(row.names(rawdata[ensembl_ids,])), res)
        
        
        # ORIGINAL iv <- match(rownames(res1),AnnotationBM$ensembl_gene_id)   # TS: deactivated gives some error

        # ORIGINAL res1$gene <- AnnotationBM[iv,'external_gene_name']    # TS: deactivated gives some error

        up <- res1$padj < 0.05 & res1$log2FoldChange > 0 #cutoffs to mess with
        flag <- rep(0, nrow(res1))
        flag[up] <- 1
        res1$up <- flag
        dn <- res1$padj < 0.05 & res1$log2FoldChange < 0 #also cutoffs to mess with
        flag1 <- rep(0, nrow(res1))
        flag1[dn] <- 1
        res1$dn <- flag1
        # ORIGINAL res1$biotype <- AnnotationBM[iv,'gene_biotype']      # TS: deactivated gives some error
        # ORIGINAL res1$description <- AnnotationBM[iv,'description']   # TS: deactivated gives some error
        res1 <- res1[order(-res1$up, -res1$dn, res1$padj),]
        colnames(res1)[1] <- "Symbol"


        out_dir = paste0(out_base, "/DE/Flu_imac/")
        if (dir.exists(out_dir)==FALSE){
          dir.create(out_dir, recursive = TRUE)
        }
        p = paste0(out_dir,i,"_pfu_",j,"_ages_",paste(k,collapse = " "),"_DE.csv")

        write.csv(res1,p,row.names = T)

      }
    } else {print(paste(i," pfu ",j," has only one age condition"))}
  }
}



savehistory(paste0(out_base, "datalog.txt"))
