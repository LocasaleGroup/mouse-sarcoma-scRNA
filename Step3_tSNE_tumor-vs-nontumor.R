rm(list=ls())
library("reticulate")
use_condaenv(condaenv="Renv", conda="/home/xiao/anaconda3/bin/conda")

library(Seurat)
library(ggplot2)

in.dir <- "QC"
out.dir <- "TSNE_cluster"
if(!dir.exists(out.dir)){
  dir.create(out.dir)
}


commArgs <- commandArgs()
sample_name <- commArgs[6]

#read the seurat obj
tumor.obj <- readRDS(file.path(in.dir,paste0(sample_name,"_filter_seurat_obj.rds")))
counts_mat <- as.matrix(tumor.obj@assays$RNA@counts)

#filter the genes
ncells <- apply(counts_mat,1,function(x) sum(x>0))
tumor.obj <- tumor.obj[ncells>10,]

## read the CNV infromation
cnv_group <- read.table(file.path(paste0(sample_name,"_infercnv_filter_out"),
                                  "infercnv.observation_groupings.txt"),
                        header=F,skip=1,row.names=1) 
tumor.obj@meta.data$cnv_group <- cnv_group[colnames(tumor.obj),"V2"]


#normalization
tumor.obj <- NormalizeData(tumor.obj)
#identification of highly variable features (feature selection)
tumor.obj<- FindVariableFeatures(tumor.obj)
top10 <- head(VariableFeatures(tumor.obj))
plot1 <- VariableFeaturePlot(tumor.obj)
plot2 <- LabelPoints(plot=plot1,points=top10,repel=TRUE)
CombinePlots(plots=list(plot1,plot2))
dev.off()

#scaling the data for PCA
all.genes <- rownames(tumor.obj)
tumor.obj<- ScaleData(tumor.obj,features = all.genes)

#PCA
tumor.obj <- RunPCA(tumor.obj, features = VariableFeatures(object = tumor.obj))
VizDimLoadings(tumor.obj,dims = 1:2, reduction="pca")
DimPlot(tumor.obj,reduction="pca")
DimHeatmap(tumor.obj, dims =1:20, cells = 500, balanced = TRUE)


#determine the dimensionality of the data
tumor.obj <- JackStraw(tumor.obj, num.replicate = 100)
tumor.obj <- ScoreJackStraw(tumor.obj, dims = 1:100)
JackStrawPlot(tumor.obj, dims =1:20)
ElbowPlot(tumor.obj,ndims=50)


#cluster the cells
tumor.obj <- FindNeighbors(tumor.obj, dims = 1:20)
tumor.obj <- FindClusters(tumor.obj, resolution = 0.5)

# tumor.obj <- RunUMAP(tumor.obj,dim=1:20)
# DimPlot(tumor.obj,reduction="umap",group.by="cnv_group")

tumor.obj <- RunTSNE(tumor.obj,dims=1:20)
p <- DimPlot(tumor.obj,reduction="tsne",group.by = "cnv_group")
ggsave(file.path(out.dir,paste0(sample_name,"_tsne.pdf")),p,width=4,height=3,device="pdf",useDingbats=FALSE)

