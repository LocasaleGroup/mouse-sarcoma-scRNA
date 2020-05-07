##comparison of tumor cells in veh and dig
rm(list=ls())
library("Seurat")
library(EnhancedVolcano)
library(matrixStats)
source("utils.R")

in.dir <- "cluster_tumor"
out.dir <- "comparison_tumor"

if(!dir.exists(out.dir)){
  dir.create(out.dir)
}

#metabolic genes
metabolism_gmt_file <- "Metabolic_genes/metabolism_kegg.gmt"
metabolism_pathways.list <- gmtPathways(metabolism_gmt_file)
metabolism_genes <- unique(unname(unlist(metabolism_pathways.list)))

#read the seurat obj
tumor.obj <- readRDS(file.path(in.dir,"tumor_obj.rds"))
tumor.obj$orig.group <- factor(substr(tumor.obj$orig.ident,1,3),levels=c("veh","dig"))



#differential test
#output the background genes
metabolism_genes <- intersect(metabolism_genes,rownames(tumor.obj))
write.table(metabolism_genes,file.path(out.dir,"background_metabolism_genes.txt"),quote=F,row.names=F,col.names=F)
write.table(rownames(tumor.obj),file.path(out.dir,"background_all_genes.txt"),quote=F,row.names=F,col.names=F)

Idents(tumor.obj) <- tumor.obj$orig.group
if(file.exists(file.path(out.dir,"diff_result.rds"))){
  diff_result <- readRDS(file.path(out.dir,"diff_result.rds"))
}else{
  diff_result <- FindMarkers(tumor.obj,logfc.threshold = -Inf,
                             min.pct = -Inf,min.cells.feature = -Inf,
                             ident.1 = "dig",ident.2="veh")
  saveRDS(diff_result,file.path(out.dir,"diff_result.rds"))
}
##calculate the cohens' d 
tumor.data <- GetAssayData(tumor.obj,assay="RNA",slot="data") %>% data.matrix %>% expm1
group_dig <- which(tumor.obj$orig.group=="dig")
group_veh <- which(tumor.obj$orig.group=="veh")
cohen_d <- cal_cohend_matrix(tumor.data,group_veh,group_dig)
diff_result$cohen_d <- cohen_d[rownames(diff_result)]

p <- EnhancedVolcano(diff_result,
                lab = rownames(diff_result),
                selectLab = F,
                x = 'cohen_d',
                y = 'p_val',
                xlab = bquote("Cohen's d"),
                ylab = bquote(~-Log[10]~italic(P~val)),
                pCutoff = 0.01,
                FCcutoff = 0.1,
                transcriptLabSize = 4.0,
                colAlpha = 0.8,
                legendVisible=F,
                title=NULL,
                subtitle=NULL,
                axisLabSize=10,
                gridlines.major = F,
                gridlines.minor = F,
                borderWidth=0.5)
#p <- p+scale_y_continuous(expand = c(0, 1))
ggsave(file.path(out.dir,"volcanoplot_all.pdf"),p,width=4.5,height=4)
#VlnPlot(tumor.obj ,"Atp5k",,pt.size=0)

#plot metabolism
p <- EnhancedVolcano(diff_result[metabolism_genes,],
                     lab = rownames(diff_result[metabolism_genes,]),
                     selectLab = F,
                     x = 'cohen_d',
                     y = 'p_val',
                     xlab = bquote("Cohen's d"),
                     ylab = bquote(~-Log[10]~italic(P~val)),
                     pCutoff = 0.01,
                     FCcutoff = 0.1,
                     transcriptLabSize = 4.0,
                     colAlpha = 0.8,
                     legendVisible=F,
                     title=NULL,
                     subtitle=NULL,
                     axisLabSize=10,
                     gridlines.major = F,
                     gridlines.minor = F,
                     borderWidth=0.5)
ggsave(file.path(out.dir,"volcanoplot_metabolism.pdf"),p,width=2.8,height=2.8,units="in")


#output the up or downregulated genes
up_result <- diff_result[which((diff_result$p_val <= 0.01)&(diff_result$cohen_d >= 0.1)),]
down_result <- diff_result[which((diff_result$p_val <= 0.01)&(diff_result$cohen_d <= -0.1)),]
write.table(rownames(up_result),file.path(out.dir,"up.genes"),quote=F,row.names = F,col.names=F)
write.table(rownames(down_result),file.path(out.dir,"down.genes"),quote=F,row.names = F,col.names=F)

##
up_metabolism <- intersect(rownames(up_result), metabolism_genes)
down_metabolism <- intersect(rownames(down_result), metabolism_genes)
write.table(up_metabolism,file.path(out.dir,"metabolism_up.genes"),quote=F,row.names = F,col.names=F)
write.table(down_metabolism,file.path(out.dir,"metabolism_down.genes"),quote=F,row.names = F,col.names=F)

up_enrich <- fisher_test(up_metabolism,metabolism_genes,metabolism_gmt_file)
down_enrich <- fisher_test(down_metabolism,metabolism_genes,metabolism_gmt_file) 
both_enrich <- fisher_test(c(up_metabolism,down_metabolism),metabolism_genes,metabolism_gmt_file)
write.csv(up_enrich,file.path(out.dir,paste0("metabolism_up_fisher.csv")))
write.csv(down_enrich,file.path(out.dir,paste0("metabolism_down_fisher.csv")))
write.csv(both_enrich,file.path(out.dir,paste0("metabolism_both_fisher.csv")))

p_up <- bubble_plot(up_enrich)
p_down <- bubble_plot((down_enrich))
p_both <- bubble_plot((both_enrich))

ggsave(file.path(out.dir,"metabolism_up_fisher.pdf"),p_up,width = 3.3,height=2.4,units="in",device="pdf",useDingbats=FALSE)
ggsave(file.path(out.dir,"metabolism_down_fisher.pdf"),p_down,width = 3.3,height=2.4,units="in",device="pdf",useDingbats=FALSE)
ggsave(file.path(out.dir,"metabolism_both_fisher.pdf"),p_both,width = 3.3,height=2.4,units="in",device="pdf",useDingbats=FALSE)

#GSEA analysis
runGSEA(tumor.obj[metabolism_genes,],test="dig",control="veh",base.name="tumor",metabolism_gmt_file,"GSEA_tumor")

##GSEA results plot
library(stringr)
pathway_name_df <- data.frame(name =names(metabolism_pathways.list),
                              stringsAsFactors = F,
                              row.names = str_to_upper(names(metabolism_pathways.list)))

file_dir <- list.files(path="GSEA_tumor",pattern = "^tumor_dig_vs_veh.Gsea.",full.names=T)
up_result_file <- list.files(path=file_dir,pattern="^gsea_report_for_dig_(.*).xls",full.names=T)
down_result_file <- list.files(path=file_dir,pattern="^gsea_report_for_veh_(.*).xls",full.names=T)
gsea_up <- read.table(up_result_file,header = T,sep="\t",row.names=1,stringsAsFactors = F)
gsea_down <- read.table(down_result_file,header = T,sep="\t",row.names=1,stringsAsFactors = F)

enrich_df <- rbind(gsea_down[1:5,], gsea_up[5:1,])
enrich_df <- enrich_df[,c(1,5,6)]
colnames(enrich_df) <- c("name","NES","pval")
enrich_df[,"name"] <- pathway_name_df[enrich_df[,"name"],]
pvals <- enrich_df[,"pval"]
pvals[which(pvals==0)] <- 0.0001
enrich_df[,"pval"] <- -log10(pvals)
enrich_df[,"pval"] <- sign(enrich_df[,"NES"]) * enrich_df[,"pval"] 
enrich_df[,"name"] <- factor(enrich_df[,"name"],levels = enrich_df[,"name"])


g <-  ggplot(enrich_df,aes(x=name,y=NES,fill=pval)) + 
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low="#0D42EC", mid="white",high="#FB291F",
                       midpoint = 0,
                       breaks=seq(-4,4,1)) +
  theme(
    axis.line = element_line(size=0.4, colour = "black"),
    # panel.grid.major = element_line(colour = "#d3d3d3",linetype = "dashed"),
    # panel.grid.minor = element_blank(),
    axis.ticks = element_line(colour = "black", size = 0.4),
    panel.border = element_blank(), panel.background = element_blank(),
    axis.text.x=element_text(colour="black", size = 8,angle=90,hjust=1,vjust=0.5),
    axis.text.y=element_text(colour="black", size = 8)) +
  theme(plot.margin = unit(rep(1,4),"lines")) 
ggsave(file.path(out.dir,"tumor_GSEA.pdf"),
       g,width = 3.3,height=5,units="in",device="pdf",useDingbats=FALSE)



