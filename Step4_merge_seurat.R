##merge the seurat objects and extract tumor and non-tumor cells
library(Seurat)

out.dir <- "merge_seurat"
if(!dir.exists(out.dir)){
  dir.create(out.dir)
}

sample_names <- c("veh1","veh3","dig1","dig2")

#read raw data
seurat_obj_list <- 
  sapply(sample_names,FUN=function(x) readRDS(file.path("QC",paste0(x,"_filter_seurat_obj.rds"))))

#merge seurat obj
merged_seurat_obj <- merge(x=seurat_obj_list[[sample_names[1]]],
                           y=seurat_obj_list[sample_names[2:length(sample_names)]],
                           add.cell.ids = sample_names,
                           project = "Digoxin")
##
print("original number of cells")
table(merged_seurat_obj$orig.ident)  #print cell numbers

#read the CNV groupings
cnv_group_list <- 
  lapply(sample_names,FUN=function(x) read.table(
    file.path(paste0(x,"_infercnv_filter_out"),"infercnv.observation_groupings.txt"),
    header=F,skip=1,row.names=1) 
    )
names(cnv_group_list) <- sample_names

#extract the tumor cells and non-tumor cells
tumor_cells <- sapply(sample_names,FUN=function(x) paste(x,rownames(subset(cnv_group_list[[x]],V2=="1")),sep="_"),USE.NAMES = F )
tumor_cells <- unlist(tumor_cells)
nontumor_cells <- sapply(sample_names,FUN=function(x) paste(x,rownames(subset(cnv_group_list[[x]],V2=="2")),sep="_"),USE.NAMES = F )
nontumor_cells <- unlist(nontumor_cells)

merged_seurat_obj_tumor <- merged_seurat_obj[,tumor_cells]
merged_seurat_obj_nontumor <- merged_seurat_obj[,nontumor_cells]

print("number of cells (tumor and non-tumor) after filtering")
table(merged_seurat_obj_tumor$orig.ident) #print tumor cell num
table(merged_seurat_obj_nontumor$orig.ident) #print nontumor cell num

#filter the genes which expressed in less then 10 cells
sels <- apply(as.matrix(merged_seurat_obj@assays$RNA@counts),1,function(x) sum(x>0)>10)
merged_seurat_obj <- merged_seurat_obj[sels,]
sels <- apply(as.matrix(merged_seurat_obj_tumor@assays$RNA@counts),1,function(x) sum(x>0)>10)
merged_seurat_obj_tumor <- merged_seurat_obj_tumor[sels,]
sels <- apply(as.matrix(merged_seurat_obj_nontumor@assays$RNA@counts),1,function(x) sum(x>0)>10)
merged_seurat_obj_nontumor <- merged_seurat_obj_nontumor[sels,]

#write to disk
saveRDS(merged_seurat_obj,file=file.path(out.dir,"merged_seurat_obj.rds"))
saveRDS(merged_seurat_obj_tumor,file=file.path(out.dir,"merged_seurat_obj_tumor.rds"))
saveRDS(merged_seurat_obj_nontumor,file=file.path(out.dir,"merged_seurat_obj_nontumor.rds"))