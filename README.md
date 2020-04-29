mouse-sarcoma-scRNA
===================

Introduction
------------
scRNA-seq analysis for studying the shifts of the tumor microenvironment upon digoxin treatment.

There are 8 main steps involved in this workflow. Four samples are involved: "veh1" and "veh3" are two replicates of control group; "dig1" and "dig2" are two replicates of digoxin treatment group.

Requirements
------------
The required R packages: Seurat, infercnv, NGCHM, ggplot2, reticulate, pheatmap, EnhancedVolcano, matrixStats, dplyr, plyr, reshape2, 

The required python package: RiboCode

Import data and quality control
-------------------------------
``` bash
Rscript Step1_QC.R
```
The gene expression profiles and annotations of cells and genes after filtering are stored as the R objects.

Copy Number Variation Analysis 
------------------------------
``` bash
python gtf_to_position_file.py 
for i in veh1 veh3 dig1 dig2
do
  Rscript Step2_inferCNV.R
done
```
The script "gtf_to_position_file.py" extracts the start and end coordinates of each gene on the genome. 
The CNV analysis will generate the graphs in Figure 4B and S3D.

Clustering of malignant and non-malignant cells 
-----------------------------------------------
``` bash
Rscript Step3_tSNE_tumor-vs-nontumor.R
```
This step will generate the graphs in Figure S3E, distinguishing the malignant cells from non-malignant cells.

Merging cells of different samples
----------------------------------
``` bash
Rscript Step4_merge_seurat.R
```
The malignant cells and non-malignant cells of all samples will be stored as the R objects separately.

Clustering the malignant cells
------------------------------
``` bash
Rscript Step5_cluster_tumor.R
```
The gene expression of malignant cells will be normalized and scaled for the downstream analysis.

Differential analysis of malignant cells
----------------------------------------
``` bash
Rscript Step6_comparison_tumor.R
```
This script will compare the gene experssion of malignant cells between control group and digoxin treatment group and perform the pathway analysis. This step will generate the graphs in Figure 5A and 5B.

Clustering of non-malignant cells
---------------------------------
``` bash
Rscript Step7_clustering_nontumor.R
```
Then non-malignant cell types will be defined in this step. This step will generate the graphs in Figure 4D.

Differential analysis of non-malignant cells
--------------------------------------------
``` bash
Rscript Step8_comparison_nontumor.R
```
This script will compare the gene experssion of non-malignant cells between control group and digoxin treatment group and perform the pathway analysis. This step will generate the graphs in Figure 5C, 5D, S5 and S6.


Contact
-------
zhengtao.xiao[at]duke.edu
