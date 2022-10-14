rm(list=ls())

library(DESeq2)

# https://figshare.com/articles/Differential_Gene_Expression/12227531
p_in = '~/Dropbox/aging_map_paper/datasets/general/resources/figshare/schaum_2020/MACA_Bulk_deseq2_results.bin'
load(p_in)


out_base = "/Users/tstoeger/Dropbox/aging_map_paper/dynamic/tstoeger/200501_schaum/"

# Prepare output folde in case it was missing
if (dir.exists(out_base)==FALSE){
  dir.create(out_base, recursive = TRUE)            
}


for (tissue in unique(ls(results_list_MACA_bulk))){
  for (comparison in unique(ls(results_list_MACA_bulk[[tissue]]))){
    to_export <- results_list_MACA_bulk[[tissue]][[comparison]]$resall
    file_name = paste0(tissue, '_', comparison, '.csv')
    p = file.path(out_base, file_name)
    write.csv(to_export,p,row.names = T)
  }
}
