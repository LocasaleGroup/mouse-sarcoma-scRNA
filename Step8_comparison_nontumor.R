#compare different condition

rm(list=ls())
library(Seurat)
library(ggplot2)
library(reshape2)
library(plyr)
library(EnhancedVolcano)
library(matrixStats)
library(stringr)
source("utils.R")

in.dir <- "cluster_nontumor_v2"
out.dir <- "comparison_nontumor"
if(!dir.exists(out.dir)){
  dir.create(out.dir)
}

#metabolism data 
metabolism_gmt_file <- "Metabolic_genes/metabolism_kegg.gmt"
metabolism_pathways.list <- gmtPathways(metabolism_gmt_file)
metabolism_genes <- unique(unname(unlist(metabolism_pathways.list)))


#read the seurat obj
nontumor.obj <- readRDS(file.path(in.dir,"nontumor.rds"))
nontumor.obj$orig.group <- factor(substr(nontumor.obj$orig.ident,1,3),levels=c("veh","dig"))
nontumor.obj$cluster_labels <- Idents(nontumor.obj)
Idents(nontumor.obj) <- nontumor.obj$orig.group

compared_cluster <- levels(nontumor.obj$cluster_labels)[1:12]


##response to digoxin for each cell type
response_result_list <- list()
###perform differential test
diff_method <- "wilcox"
response_logFC_df <- as.data.frame(matrix(NA,nrow=nrow(nontumor.obj),
                                          ncol=length(compared_cluster),
                                          dimnames = list(rownames(nontumor.obj),compared_cluster)))
response_pval_adj_df <- response_logFC_df
response_pval_adj_df_d <- response_logFC_df
response_pval_df <- response_logFC_df
response_pval_df_d <- response_logFC_df
response_corhen_d <- response_logFC_df

#if the file exists
if(file.exists(file.path(out.dir,"diff_result_list_wilcox.rds"))){
  diff_result_list <- readRDS(file.path(out.dir,paste0("diff_result_list_",diff_method,".rds")))
}else{
  diff_result_list <- list()
}

for(i in compared_cluster){
  cluster_i_obj <- nontumor.obj[,nontumor.obj$cluster_labels==i]
  if(is.null(diff_result_list[[i]])){
    diff_result <- FindMarkers(cluster_i_obj,
                               test.use=diff_method,
                               logfc.threshold = -Inf,
                               min.pct = -Inf,min.cells.feature = -Inf,
                               ident.1 = "dig",ident.2="veh")
    diff_result_list[[i]] <- diff_result
  }else{
    diff_result <- diff_result_list[[i]]
  }
  
  response_logFC_df[rownames(diff_result),i] <- diff_result$avg_logFC
  pval <- diff_result$p_val
  uniq_pval <- unique(pval[!is.na(pval)])
  min_pval <- sort(uniq_pval[which(uniq_pval!=0)])[1]
  pval[pval==0] <- min_pval
  #replace the zero as the minimal value
  response_pval_df[rownames(diff_result),i] <- pval
  response_pval_df_d[rownames(diff_result),i] <- -log10(pval) * sign(diff_result$avg_logFC)
  response_pval_adj_df[rownames(diff_result),i] <- diff_result$p_val_adj
  response_pval_adj_df_d[rownames(diff_result),i] <- -log10(diff_result$p_val_adj) * sign(diff_result$avg_logFC)
  
  #calculate the Cohen's d
  cluster_i_data <- GetAssayData(cluster_i_obj,assay="RNA",slot="data") %>% data.matrix 
  group_dig <- which(cluster_i_obj$orig.group=="dig")
  group_veh <- which(cluster_i_obj$orig.group=="veh")
  cohen_d <- cal_cohend_matrix(cluster_i_data,group_veh,group_dig)
  response_corhen_d[rownames(cluster_i_data),i] <- cohen_d
}
response_result_list <- list(response_logFC_df,response_pval_adj_df,response_pval_adj_df_d,
                                  response_pval_df,response_pval_df_d,response_corhen_d)
if(!file.exists(file.path(out.dir,"diff_result_list_wilcox.rds"))){
  saveRDS(diff_result_list,file.path(out.dir,"diff_result_list_wilcox.rds"))
}


#########################Compare the TME in Veh and Dig groups
choice <- "metabolic_genes" #or "all_genes"
if(choice == "metabolic_genes"){
  select_genes <-intersect(metabolism_genes,rownames(nontumor.obj))
  select_name <- "metabolic"
}else{
  select_genes <- rownames(nontumor.obj)
  select_name <- "all"
}
#####################################################################################
#output the background genes
write.table(select_genes,file.path(out.dir,paste0("background_",select_name,"_genes.txt")),quote=F,row.names=F,col.names=F)

####plot and enrichment analysis
response_logFC_df <- response_result_list[[1]]
response_pval_adj_df <- response_result_list[[2]]
response_pval_adj_df_d <- response_result_list[[3]]
response_pval_df <- response_result_list[[4]]
response_pval_df_d <- response_result_list[[5]]
response_corhen_d <- response_result_list[[6]]
#plot
#plot the genes: cohens's d >0.1 and pvalue <0.05 and do differential test
response_data <- response_corhen_d
response_data <- response_data[select_genes,]
figname <- paste0(select_name,"_wilcox_cohen_sig.pdf") #_logFC, _corhen.pdf , _pval.pdf, _padj.pdf
####plot differential genes~~~~~
plot_data <- melt(response_data)
pval_data <- melt(response_pval_df[select_genes,])
plot_data$color <- 1
plot_data[which((plot_data$value >= 0.1) & (pval_data$value <= 0.05)),"color"] = 2 #
plot_data[which((plot_data$value <= -0.1) & (pval_data$value <= 0.05)),"color"] = 3 
plot_data$color <- factor(plot_data$color,levels=c(2,1,3))

g <- ggplot(plot_data,aes(x=variable,y=value,order=color)) +
  labs(y="Cohens' d",x=NULL)+
  geom_jitter(aes(colour=color),na.rm=TRUE,size=0.5,alpha=0.7) +
  #stat_summary(fun.y = median,geom="point",size=1,color="red")+
  scale_color_manual(values=c("firebrick3", "#ede1cc","dodgerblue4")) +
  #geom_hline(yintercept = 0,col="grey") +
  #scale_y_continuous(breaks = seq(-1.8,1,0.2)) +
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 8,angle=45,hjust=1,vjust=1),
        axis.text.y=element_text(colour="black", size = 8),
        axis.line=element_line(size=0.45,color="black"),
        axis.ticks = element_line(colour = "black",size=0.45),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"))
ggsave(file.path(out.dir,figname),g,height=2.5,width=3,useDingbats=FALSE)
####plot end~~~~~

#pathway analysis for differential genes
group_enrich_results <- list()
for(i in compared_cluster){
  #output the differential genes
  up_result <- response_pval_df[which((response_pval_df[,i] <= 0.05)&(response_corhen_d[,i] >= 0.1)),]
  down_result <- response_pval_df[which((response_pval_df[,i] <= 0.05)&(response_corhen_d[,i] <= -0.1)),]
  if(select_name == "metabolic"){
    up_result <- up_result[intersect(select_genes,rownames(up_result)),]
    down_result <- down_result[intersect(select_genes,rownames(down_result)),]
  }
  write.table(rownames(up_result),file.path(out.dir,paste0(select_name,"_",diff_method,"_cluster_",i,"_up.genes")),quote=F,row.names = F,col.names=F)
  write.table(rownames(down_result),file.path(out.dir,paste0(select_name,"_",diff_method,"_cluster_",i,"_down.genes")),quote=F,row.names = F,col.names=F)
  
  
  if(select_name == "metabolic"){
    ##using fisher exact test to perform enrichment test
    up_enrich <- fisher_test(rownames(up_result),metabolism_genes,metabolism_gmt_file)
    down_enrich <- fisher_test(rownames(down_result),metabolism_genes,metabolism_gmt_file) 
    both_enrich <- fisher_test(c(rownames(up_result),rownames(down_result)),
                               metabolism_genes,metabolism_gmt_file)
    write.csv(up_enrich,file.path(out.dir,paste0(select_name,"_",diff_method,"_cluster_",i,"_up_fisher.csv")))
    write.csv(down_enrich,file.path(out.dir,paste0(select_name,"_",diff_method,"_cluster_",i,"_down_fisher.csv")))
    write.csv(both_enrich,file.path(out.dir,paste0(select_name,"_",diff_method,"_cluster_",i,"_both_fisher.csv")))
    group_enrich_results[[i]] <- both_enrich    
  }else{
    print("please do KEGG analysis on metascape website")
  }  
}

##keep top 3 pathways
group_enrich_results <- lapply(group_enrich_results, FUN=function(x) x[order(x$pval),][1:min(5,nrow(x)),])
##add the cell type and pathway name
group_enrich_results <- lapply(names(group_enrich_results), 
                               FUN=function(x) cbind(rownames(group_enrich_results[[x]]),group_enrich_results[[x]],rep(x,nrow(group_enrich_results[[x]]))) )
#rbind
group_enrich_results <- do.call(rbind,group_enrich_results) 
colnames(group_enrich_results) <- c("pathway","pval","hit","rich","celltype")

####!!!! replace the pval that smaller than 0.001 as the 0.001
pvals <- group_enrich_results$pval
pvals[which(pvals<0.001)] <- 0.001
group_enrich_results$pval <- pvals
p <- group_bubble_plot(group_enrich_results)
ggsave(file.path(out.dir,"pathway_enrich.pdf"),
       p,width=9.8,height=5,units="in",device="pdf",useDingbats=F)


#####
#GSEA analysis
if(select_name == "metabolic"){
  for(i in compared_cluster){
    cluster_i_obj <- nontumor.obj[select_genes,nontumor.obj$cluster_labels==i]
    runGSEA(cluster_i_obj,test="dig",control="veh",base.name=i,metabolism_gmt_file,"GSEA_notumor")
  }
}

##GSEA results plot
pathway_name_df <- data.frame(name =names(metabolism_pathways.list),
                              stringsAsFactors = F,
                              row.names = str_to_upper(names(metabolism_pathways.list)))
for(i in compared_cluster){
  file_dir <- list.files(path="GSEA_notumor",pattern = paste0("^",i,"_dig_vs_veh.Gsea."),full.names=T)
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
  ggsave(file.path(out.dir,paste0(i,"_GSEA.pdf")),
         g,width = 3.3,height=5,units="in",device="pdf",useDingbats=FALSE)
}

