#install some dependence packages
install.packages("Seurat")#‘5.0.3’
install.packages("ggpubr")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicRanges")

devtools::install_github("sulab-wmu/scPagwas")
library(ggpubr)
library(Seurat)
library(GenomicRanges)
library(scPagwas)
library(Seurat)
library(scPagwas)
cell_path="G:/07.STAD_subtype/analysis0507_v7/07.scRNA-seq/harmony_0.5_30/markers_by_self"
Single_data<-readRDS(file.path(cell_path,"Seurat_Anno_cells_by_myself.rds"))
metadata=Single_data@meta.data
Idents(Single_data) <- Single_data$cell_type
setwd("G:/07.STAD_subtype/analysis0507_v7/26.scPagwas")
load("mp_immune_genelst.RData")
Pagwas_data<-scPagwas_main2(Pagwas = NULL,
                            gwas_data ="finngen_r7_c3_stomach_exallc.txt", # The GWAS Summary statistics files 
                            Single_data =Single_data,
                            output.prefix="stomach", # the prefix name for output files
                            output.dirs="scPagwas_sumu_immnue",# the directory file's name for output
                            block_annotation = block_annotation,
                            assay="SCT", # the assays for scRNA-seq data to use.
                            Pathway_list=mp_immune_genelst,# pathway list is provided by package, including gene symbols.
                            n.cores=1,
                            iters_singlecell = 100,
                            chrom_ld = chrom_ld,# The LD data is provided by package.
                            celltype=T# Whether to run the celltype process.
)
saveRDS(Pagwas_data,"Pagwas_data.rds")
##--结果可视化---
rm(list=ls())
setwd("G:/07.STAD_subtype/analysis0507_v7/26.scPagwas")
source("scpagwas_function.R")
Pagwas_data=readRDS("Pagwas_data.rds")
setwd("G:/07.STAD_subtype/analysis0507_v7/26.scPagwas/scPagwas_sumu_immnue")
dir.create("figure")
library(ggplot2)
##---细胞与性状相关的显著性，p值是合并的p值----
cell_p=read.csv("stomach_Merged_celltype_pvalue.csv",row.names = 1)
# 添加 -log10(p-value)
cell_p$log_p_value <- -log10(cell_p$pvalue)
# 绘制棒棒图
p1 <- ggplot(cell_p, aes(x = reorder(celltype, log_p_value), y = log_p_value)) +
  geom_bar(stat = "identity", fill = "skyblue",width=0.6) +
  labs(x = "", y = "-log10(p-value)", title = "Trait Related Cell Significance") +
  theme_minimal() +
  coord_flip()  # 水平显示，便于阅读
pdf("cell_type_p.pdf",width=3,height=2)
p1
dev.off()
#The merged singlecell pvalue to celltype pvalue

##---重要得分的可视化-----
scPagwas_Visualization1(Single_data=Pagwas_data,
                       p_thre = 0.05,
                       FigureType = "tsne",
                       width = 7,
                       height = 7,
                       lowColor = "white", 
                       highColor = "red",
                       output.dirs="figure",
                       size = 0.5,
                       do_plot = FALSE)
##---重要基因的可视化----
CpG_path="G:/07.STAD_subtype/analysis0507_v7/04.OS/CpG"
meth_path="G:/07.STAD_subtype/analysis0507_v7/02.meth_diff"
CpG=read.table(file.path(CpG_path,"single_cox_sig.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
CpG_gene=read.table(file.path(meth_path,"sig_cg_gene_point.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
CpG_gene=CpG_gene[CpG_gene$cg %in% CpG$Characteristics,]
TF_path="G:/07.STAD_subtype/analysis0507_v7/04.OS/TF"
TF=read.table(file.path(TF_path,"single_cox_sig.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
TF_sig=unique(TF$Characteristics)
gene_heri_cor[which(gene_heri_cor$PCC>0.1 & gene_heri_cor$pvalue<0.05),]
heritability_cor_scatterplot1(gene_heri_cor=Pagwas_data@misc$PCC,
                             topn_genes_label=10,
                             important=TF_sig,
                             color_low="#035397",
                             color_high ="#F32424",
                             color_mid = "white",
                             text_size=2,
                             do_plot=T,
                             max.overlaps =20,
                             width =6,
                             height = 6,
                             figurenames="TF_PCC.pdf")
CpG_gene <- CpG_gene$gene
heritability_cor_scatterplot1(gene_heri_cor=Pagwas_data@misc$PCC,
                              topn_genes_label=10,
                              important=CpG_gene,
                              color_low="#035397",
                              color_high ="#F32424",
                              color_mid = "white",
                              text_size=2,
                              do_plot=T,
                              max.overlaps =20,
                              width =6,
                              height = 6,
                              figurenames="CpG_gene_PCC.pdf")


