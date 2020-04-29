library(Seurat)
library(infercnv)
library(NGCHM)
library("stringr")

args <- commandArgs()
sample_name <- args[6] # "veh1","veh3","dig1","dig2" 
min_cells <- 10
min_features <- 1800
out.dir <- "infercnv_input"
if(!dir.exists(out.dir)){
  dir.create(out.dir)
}

print(paste("starting",sample_name))

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = file.path(out.dir,paste0(sample_name,".merged.counts.matrix")),
                                     annotations_file = file.path(out.dir,paste0(sample_name,".merged.cell_annot.txt")),
                                     delim="\t",
                                     gene_order_file = "mm10_gene_pos.txt",
                                     ref_group_names = c("muscle14","muscle15"),
                                     chr_exclude = c("X","Y","M"))

infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir=paste0(sample_name,"_infercnv_filter_out"),
                              denoise=T,
                              hclust_method="ward.D2",
                              k_obs_groups=2,
                              HMM=F,
                              num_threads = 16)
