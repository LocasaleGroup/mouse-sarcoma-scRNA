rm(list=ls())
library("reticulate")
use_condaenv(condaenv="Renv", conda="/home/xiao/anaconda3/bin/conda")

library(Seurat)
library(ggplot2)
library(pheatmap)
library(stringr)
library(dplyr) 
library(reshape2)
source("utils.R")


in.dir <- "merge_seurat"
out.dir <- "cluster_nontumor_v2"
if(!dir.exists(out.dir)){
  dir.create(out.dir)
}
sample_names <- c("veh1","veh3","dig1","dig2")
sample_group <- factor(c("veh","veh","dig","dig"),levels=c("veh","dig"))

#read the seurat obj
nontumor.obj <- readRDS(file.path(in.dir,"merged_seurat_obj_nontumor.rds"))
nontumor.obj$orig.group <- factor(substr(nontumor.obj$orig.ident,1,3),levels=c("veh","dig"))

#normalization
nontumor.obj <- NormalizeData(nontumor.obj)
#identification of highly variable features (feature selection)
nontumor.obj<- FindVariableFeatures(nontumor.obj,nfeatures=2000)
top10 <- head(VariableFeatures(nontumor.obj))
plot1 <- VariableFeaturePlot(nontumor.obj)
plot2 <- LabelPoints(plot=plot1,points=top10,repel=TRUE)
CombinePlots(plots=list(plot1,plot2))
dev.off()

#scaling the data for PCA
all.genes <- rownames(nontumor.obj)
nontumor.obj<- ScaleData(nontumor.obj,features = all.genes)

#PCA
nontumor.obj <- RunPCA(nontumor.obj, features = VariableFeatures(object = nontumor.obj),npcs=100)
VizDimLoadings(nontumor.obj,dims = 1:5,reduction="pca")
DimPlot(nontumor.obj,reduction="pca")
DimHeatmap(nontumor.obj, dims =1:20, cells = 500, balanced = TRUE)

#determine the dimensionality of the data
nontumor.obj <- JackStraw(nontumor.obj,dims=100)
nontumor.obj <- ScoreJackStraw(nontumor.obj, dims = 1:100)
JackStrawPlot(nontumor.obj, dims =1:100)
ElbowPlot(nontumor.obj,ndims=100)

select_pcs <- 85

#cluster the cells
nontumor.obj <- FindNeighbors(nontumor.obj, dims = 1:select_pcs)
nontumor.obj <- FindClusters(nontumor.obj,resolution = 0.15)

perplexity <- max(30,floor(ncol(nontumor.obj)/100))
eta <- max(200,floor(ncol(nontumor.obj)/12))

nontumor.obj <- RunTSNE(nontumor.obj,dims=1:select_pcs,tsne.method="Rtsne",perplexity=perplexity,eta=eta)
p <- DimPlot(nontumor.obj,reduction="tsne")
ggsave(file.path(out.dir,"nontumor_tsne.pdf"),p,width=5,height=4,device="pdf",useDingbats=FALSE)


### assign the cell types
new.cluster.ids <- c("M2", "M1", "M2#", "T_cell", "NK", 
                     "Epithelial","Osteoclasts","Fibroblast",
                     "Microglia","Neutrophil","DC","B_cell","Unknown")
names(new.cluster.ids) <- levels(nontumor.obj)
nontumor.obj <- RenameIdents(nontumor.obj, new.cluster.ids)



###output the seurat obj 
saveRDS(nontumor.obj,file.path(out.dir,"nontumor.rds"))
#saveRDS(nontumor.markers,file.path(out.dir,"nontumor.markers"))
