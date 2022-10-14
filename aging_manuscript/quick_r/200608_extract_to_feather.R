library(Seurat)
library(feather)

store_as_feather <- function(data_to_export, p_out) {
  # p_out = file.path(
  #   "/Users/tstoeger/Desktop", 
  #   paste(c(outputname, ".feather"), collapse = "")
  # )
  df <- as.data.frame(as.matrix(data_to_export))
  df <- cbind(feather_index = rownames(df), df)  # feather doesn't save row indices, which could lead to loss of identifier
  write_feather(df, p_out)
}

out_main_folder = '/Users/tstoeger/Dropbox/aging_map_paper/dynamic/tstoeger/200608_single_cell'
dir.create(out_main_folder, showWarnings = FALSE)

sub_name = 'tabula_muris_10x'
in_folder = '/Users/tstoeger/Dropbox/aging_map_paper/datasets/general/resources/figshare/pisco_2019_8273102'
organs = c(
  'Bladder_droplet',
  'Fat_droplet',
  'Heart_droplet',
  'Kidney_droplet',
  'Large_Intestine_droplet',
  'Limb_Muscle_droplet',
  'Liver_droplet',
  'Lung_droplet',
  'Mammary_Gland_droplet',
  'Pancreas_droplet',
  'Skin_droplet',
  'Spleen_droplet',
  'Thymus_droplet',
  'Tongue_droplet',
  'Trachea_droplet'
)


out_folder = file.path(out_main_folder, sub_name)
dir.create(out_folder, showWarnings = FALSE)

for (organ in organs){
  p_in <- file.path(in_folder, paste0(organ, '.h5ad'))
  seurat_object <- ReadH5AD(p_in)
  transcript_counts <- seurat_object@assays$RNA@counts
  p_out <- file.path(out_folder, paste0(organ, '.feather'))
  store_as_feather(transcript_counts, p_out)
}


sub_name = 'tabula_muris_facs'
in_folder = '/Users/tstoeger/Dropbox/aging_map_paper/datasets/general/resources/figshare/pisco_2019_8273102'
organs = c(
  'Bladder_facs',
  'Brain_Myeloid_facs',
  'Brain_Non-Myeloid_facs',
  'Diaphragm_facs',
  'Heart_facs',
  'Kidney_facs',
  'Large_Intestine_facs',
  'Limb_Muscle_facs',
  'Liver_facs',
  'Lung_facs',
  'Mammary_Gland_facs',
  'Marrow_facs',
  'Pancreas_facs',
  'Skin_facs',
  'Spleen_facs',
  'Thymus_facs',
  'Tongue_facs',
  'Trachea_facs'
)


out_folder = file.path(out_main_folder, sub_name)
dir.create(out_folder, showWarnings = FALSE)

for (organ in organs){
  p_in <- file.path(in_folder, paste0(organ, '.h5ad'))
  seurat_object <- ReadH5AD(p_in)
  transcript_counts <- seurat_object@assays$RNA@counts
  p_out <- file.path(out_folder, paste0(organ, '.feather'))
  store_as_feather(transcript_counts, p_out)
}


sub_name = 'calico'
in_folder = '/Users/tstoeger//Dropbox/aging_map_paper/datasets/general/resources/publications/kimmel_2019/download_from_calico_190609/'
organs = c(
  'lung', 'kidney', 'spleen'
)


out_folder = file.path(out_main_folder, sub_name)
dir.create(out_folder, showWarnings = FALSE)

for (organ in organs){
  p_in <- file.path(in_folder, paste0(organ, '.h5ad'))
  seurat_object <- ReadH5AD(p_in)
  transcript_counts <- seurat_object@assays$RNA@counts
  p_out <- file.path(out_folder, paste0(organ, '.feather'))
  store_as_feather(transcript_counts, p_out)
}



