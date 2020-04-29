library(Seurat)
library(ggplot2)
library("stringr")

sample_names <- c("veh1","veh3","dig1","dig2")
min_gene_num <- 1800
min_cell <- 10
out.dir <- "QC"
if(!dir.exists(out.dir)){
  dir.create(out.dir)
}


#############
seurat.obj.list <- list()
for(sample_name in sample_names){
  raw_data <- Read10X(data.dir = file.path("../filtered_feature_bc_matrix/",sample_name))
  seurat.obj <- CreateSeuratObject(counts = raw_data,project=sample_name,min.cells = 10)
  saveRDS(seurat.obj,file.path(out.dir,paste0(sample_name,"_raw_seurat_obj.rds")))
  seurat.obj.list[[sample_name]] <- seurat.obj
}

#QC: 
#1. number of detected genes
number_of_genes <- sapply(sample_names,FUN=function(x) cbind(rep(x,ncol(seurat.obj.list[[x]])),seurat.obj.list[[x]]@meta.data$nFeature_RNA))
number_of_genes_df <- do.call(rbind,number_of_genes)
number_of_genes_df <- data.frame(number_of_genes_df,stringsAsFactors =F)
colnames(number_of_genes_df) = c("sample","num")
number_of_genes_df$sample <- factor(number_of_genes_df$sample,levels=sample_names)
number_of_genes_df$num <- as.numeric(number_of_genes_df$num)
p <- ggplot(number_of_genes_df,aes(sample,num)) +
  geom_violin(trim = T, size=0.5,show.legend = F,width=1.0,fill="#3282bd",color="#3282bd")+
  labs(title=NULL,x=NULL, y = "Number of genes")+
  geom_boxplot(width=0.1,outlier.size=0.6)+
  geom_hline(yintercept = min_gene_num,color="blue", linetype="dashed", size=0.5) +
  scale_y_continuous(breaks = round(seq(0, max(number_of_genes_df$num), by = 1000)),expand = c(0, 0)) +
  expand_limits(y = 0) +
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 8),
        axis.text.y=element_text(colour="black", size = 8),
        axis.line=element_line(size=0.2,color="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm")) 
ggsave(file.path(out.dir,"raw_gene_num_violin.pdf"),p,width=4,height=2,units="in",device="pdf",useDingbats=FALSE)


p <- ggplot(number_of_genes_df,aes(x=num)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white",binwidth=250)+
  geom_density(alpha=.2, fill="#FF6666") +
  scale_x_continuous(breaks = round(seq(0, max(number_of_genes_df$num), by = 1000))) +
  geom_vline(xintercept = min_gene_num,color="blue", linetype="dashed", size=0.5) +
  facet_wrap(~sample,ncol=2) +
  theme_minimal() +
  theme(axis.text.x=element_text(colour="black", size = 8,hjust=1,angle = 90),
        axis.text.y=element_text(colour="black", size = 8),
        panel.border = element_blank(), panel.background = element_blank()) 
ggsave(file.path(out.dir,"raw_gene_num_hist.pdf"),p,width=4,height=3.5,units="in",device="pdf",useDingbats=FALSE)

#Filter and save
seurat.obj.filter.list <- sapply(sample_names,FUN=function(x) subset(seurat.obj.list[[x]], nFeature_RNA > min_gene_num))
lapply(sample_names,FUN=function(x) saveRDS(seurat.obj.filter.list[[x]],
                                            file=file.path(out.dir,paste0(x,"_filter_seurat_obj.rds"))
                                            )
       )
