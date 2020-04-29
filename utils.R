gmtPathways <- function(gmt.file) {
  pathwayLines <- strsplit(readLines(gmt.file), "\t")
  pathways <- lapply(pathwayLines, tail, -2)
  names(pathways) <- sapply(pathwayLines, head, 1)
  pathways
}

##calculate Cohen's d
#https://toptipbio.com/cohens-d/
#1. if two groups have same sample size: (m2-m1) / sqrt((sd1^2 + sd2^2) / 2)
#2. if two groups have different sample size: (m2 - m1) / sqrt( ((n1-1)*(sd1^2) + (n2-1)*(sd2^2)) / (n1+n2-2) )
cal_cohend <- function(values1, values2){
  n1 <- length(values1)
  n2 <- length(values2)
  m1 <- mean(values1)
  m2 <- mean(values2)
  sd1 <- sd(values1)
  sd2 <- sd(values2)
  
  if(n1 == n2){
    pooled_sd <- sqrt( (sd1^2 + sd2^2) / 2 ) 
  }else{
    pooled_sd <- sqrt( ( (n1-1)*(sd1^2) + (n2-1)*(sd2^2) ) / (n1+n2-2) )
  }
  
  cohen_d <- (m2 - m1) / pooled_sd
  
  return(cohen_d)
}

#calculate cohen's d from matrix data
cal_cohend_matrix <- function(values_matrix, group1, group2){
  n1 <- length(group1)
  n2 <- length(group2)
  
  if(any(c(n1,n2) == 0)){
    stop("error,either group1 or group2 has no samples")
  }
  
  group1_means <- rowMeans(values_matrix[,group1])
  group2_means <- rowMeans(values_matrix[,group2])
  
  group1_sds <- matrixStats::rowSds(values_matrix[,group1])
  group2_sds <- matrixStats::rowSds(values_matrix[,group2])
  
  if(n1 == n2){
    pooled_sds <- sqrt( (group1_sds + group2_sds) / 2 )
  }else{
    pooled_sds <- sqrt( ((n1-1)*(group1_sds^2) + (n2-1)*(group2_sds^2)) / (n1+n2-2) )
  }
  
  cohen_ds <- (group2_means - group1_means) / pooled_sds
  cohen_ds
} 


###fisher exact test,enrichment analysis
fisher_test <- function(input,background,gmtfile){
  input <- unique(input)
  background <- unique(background)
  numBackgroud <- length(background)
  numInput <- length(input)
  pathway.list <- gmtPathways(gmtfile)
  res <- data.frame(pval = rep(1,length(pathway.list)),
                    hit  = rep(1,length(pathway.list)),
                    rich = rep(1,length(pathway.list)),
                    row.names = names(pathway.list))
  
  for(p in names(pathway.list)){
    pathway_background <- intersect(background,pathway.list[[p]])
    numM <- length(pathway_background) 
    numHits <- length(intersect(input,pathway_background))
    
    contMat <- cbind(sig = c(numHits, numM - numHits), 
                     notSig = c(numInput-numHits, numBackgroud - numM - numInput + numHits))
    row.names(contMat) <- c("anno", "notAnno")
    if(all(contMat ==0)){
      p.value <- 1
    }else{
      p.value <- fisher.test(contMat,alternative="greater")$p.value
    }
    rich_factor <- numHits / numM
    res[p,"pval"] <- p.value
    res[p,"hit"] <- numHits
    res[p,"rich"] <- rich_factor
  }
  
  #filter: min Hits: 3, min Pvalue: 0.1
  select_hit <- res$hit >=3
  select_pval <- res$pval < 0.1
  res <- res[(select_hit & select_pval),]
  return(res[order(res$pval),,drop=F])
}


#bubble plot for all cell types
group_bubble_plot <- function(plot_data){
  library(ggplot2)
  
  library(RColorBrewer)   # for brewer.pal(...)
  plot_data$pval <- -log10(plot_data$pval)
  
  p <- ggplot(plot_data, aes(x = rich, y = pathway, size = hit, fill = pval)) +
    geom_point(shape=21) +
    expand_limits(x = 0, y = 0)+
    coord_cartesian(clip = 'off')+
    labs(size = "i", fill = "p",
         x = "Rich factor", y = "Metabolic pathways") +
    scale_size(range = c(1,6)) +
    scale_fill_gradient(low = "lightblue", high = "red",breaks= -log10(c(0.05,0.01,0.001)))  +
    facet_wrap(~celltype,nrow=1) +
    theme(
      axis.line = element_line(size=0.4, colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(colour = "black", size = 0.4),
      panel.border = element_blank(), panel.background = element_blank(),
      strip.text.x = element_text(size = 5),
      axis.text.x=element_text(colour="black", size = 6, angle=90,hjust=1,vjust=0.5),
      axis.text.y=element_text(colour="black", size = 6)) +
    theme(plot.margin = unit(rep(1,4),"lines"))
  return(p)
  
  
}

#########################################################################################################################
#GSEA
runGSEA<-function(seurat_obj,test="dig",control="veh",base.name="t",gmt.file,outdir){
  testname<-paste(base.name,test,'vs',control,sep='_')
  testname<-gsub(' ','_',testname)
  
  # save expression matrix
  scaled_data <- GetAssayData(seurat_obj,assay="RNA",slot="data")
  write.table(rbind(c('symbols',colnames(seurat_obj)),
                    cbind(rownames(seurat_obj),as.matrix(scaled_data))),
              file='expr.txt',
              quote=F,
              sep='\t',
              col.names=F,
              row.names=F)
  
  #save the cls file
  pheno<-as.character(seurat_obj$orig.group)
  con<-file('pheno.cls',open='w')
  write(paste(length(pheno),'2 1'),con)
  write(paste('# ',test, ' ',control,sep=''),con)
  classes<-''
  for (i in 1:length(pheno)){
    classes<-paste(classes,pheno[i])
  }
  write(classes,con)
  close(con)         
  
  #call java gsea version 
  command <- paste('java -Xmx512m -cp gsea-3.0.jar xtools.gsea.Gsea -res expr.txt -cls pheno.cls#',test,'_versus_',control,' -gmx ',gmt.file,
                   ' -collapse false -nperm 10000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label ',testname, 
                   ' -metric Diff_of_Classes -sort real -order descending -include_only_symbols false -make_sets true -median false -num 100',
                   ' -plot_top_x 20 -rnd_seed 123456 -save_rnd_lists false -set_max 10000 -set_min 5 -zip_report false -out ', outdir, ' -gui false',sep='')
  if(get_os() == "win"){
    system(command,show.output.on.console=F)
    #print(command)
  }else{
    system(command)
    #print(command)
  }
  
  #unlink(c('pheno.cls','expr.txt'))
}

## get_os is copied from https://github.com/r-lib/rappdirs/blob/master/R/utils.r#L1
get_os <- function() {
  if (.Platform$OS.type == "windows") { 
    "win"
  } else if (Sys.info()["sysname"] == "Darwin") {
    "mac" 
  } else if (.Platform$OS.type == "unix") { 
    "unix"
  } else {
    stop("Unknown OS")
  }
}



