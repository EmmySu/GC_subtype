##---画每个数据集的特征的热图---------
rm(list=ls())
library(pheatmap)
##--function---#
pheatmap_clincial <- function(score,network1,cluster){
  cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
  cluster=cluster[order(cluster$cluster),]
  score=score[,cluster$sample]
  score=merge(network1,score,by.x="TF",by.y="row.names")
  rownames(score)=score$regulon
  score=score[,-c(1,2)]
  clustercolor=c("#61b3de","#ffee93") 
  names(clustercolor)=c("IM","EP") #类型颜色
  ann_colors=list(cluster=clustercolor)
  annotation_col=data.frame(cluster=cluster$cluster)
  row.names(annotation_col)=cluster$sample
  p=pheatmap(as.matrix(score),cluster_rows = T,cluster_cols = T,
             color=colorRampPalette(c("navy","white","firebrick3"))(100),
             fontsize_number =6,
             show_colnames = F,border_color = NA,scale = "row",show_rownames =T,
             annotation_col = annotation_col,# annotation_row = annotation_row,
             annotation_colors = ann_colors)
  p1=pheatmap(score,cluster_rows =T,cluster_cols = F,
              color=colorRampPalette(c("navy","white","firebrick3"))(100),
              show_colnames = F,border_color = NA,scale = "row",show_rownames =T,
              fontsize_number =6,
              annotation_col = annotation_col,# annotation_row = annotation_row,
              annotation_colors = ann_colors)
  p_list <- list(p,p1)
  return(p_list)
}
###--read data---#
setwd("D:/project/07.STAD_subtype/analysis0507_v7/17.cluster_valiation")
GSE211704=read.table("GSE211704/ConsensusClusterPlus10-100-km-euclidean/consensusClass2.txt",sep="\t",header=F,stringsAsFactors = F,quote = "")
colnames(GSE211704)=c("sample","cluster")
GSE211704$cluster=paste0("cluster",GSE211704$cluster)
GSE211704$cluster[GSE211704$cluster=="cluster1"]="IM"
GSE211704$cluster[GSE211704$cluster=="cluster2"]="EP"
GSE211704_score=read.table("GSE211704/GSE211704_score.xls",sep="\t",header=T,stringsAsFactors = F,quote="")

GSE127857=read.table("GSE127857/ConsensusClusterPlus10-100-km-euclidean/consensusClass2.txt",sep="\t",header=F,stringsAsFactors = F,quote = "")
colnames(GSE127857)=c("sample","cluster")
GSE127857$cluster=paste0("cluster",GSE127857$cluster)
GSE127857$cluster[GSE127857$cluster=="cluster2"]="IM"
GSE127857$cluster[GSE127857$cluster=="cluster1"]="EP"
GSE127857_score=read.table("GSE127857/GSE127857_tumor_score.xls",sep="\t",header=T,stringsAsFactors = F,quote="")

GSE207846=read.table("GSE207846/ConsensusClusterPlus10-100-km-euclidean/consensusClass2.txt",sep="\t",header=F,stringsAsFactors = F,quote = "")
colnames(GSE207846)=c("sample","cluster")
GSE207846$cluster=paste0("cluster",GSE207846$cluster)
GSE207846$cluster[GSE207846$cluster=="cluster2"]="IM"
GSE207846$cluster[GSE207846$cluster=="cluster1"]="EP"
GSE207846_score=read.table("GSE207846/GSE207846_score.xls",sep="\t",header=T,stringsAsFactors = F,quote="")

GSE164988=read.table("GSE164988/ConsensusClusterPlus10-100-km-euclidean/consensusClass2.txt",sep="\t",header=F,stringsAsFactors = F,quote = "")
colnames(GSE164988)=c("sample","cluster")
GSE164988$cluster=paste0("cluster",GSE164988$cluster)
GSE164988$cluster[GSE164988$cluster=="cluster1"]="IM"
GSE164988$cluster[GSE164988$cluster=="cluster2"]="EP"
GSE164988_score=read.table("GSE164988/GSE164988_tumor_score.xls",sep="\t",header=T,stringsAsFactors = F,quote="")

combat=read.table("combat/ConsensusClusterPlus10-100-km-euclidean/consensusClass2.txt",sep="\t",header=F,stringsAsFactors = F,quote = "")
colnames(combat)=c("sample","cluster")
combat$cluster=paste0("cluster",combat$cluster)
combat$cluster[combat$cluster=="cluster1"]="IM"
combat$cluster[combat$cluster=="cluster2"]="EP"
combat_score=read.table("combat/combat_score.xls",sep="\t",header=T,stringsAsFactors = F,quote="")

network_path="D:/project/07.STAD_subtype/analysis0507_v4/05.cox_TF_cox_CpG_cor/network_regulon"
network=read.table(file.path(network_path,"TF_sig_cox_regulon_cor.txt"),sep='\t',header=T,stringsAsFactors = F,quote="")
network=network[order(network$number,network$abs,decreasing = T),]
network1=network[,c("TF","regulon")]

p_list_GSE211704=pheatmap_clincial(GSE211704_score,network1,GSE211704)
pdf("GSE211704/GSE211704_pheatmap.pdf",width=5,height=5)
print(p_list_GSE211704[[1]])
dev.off()
pdf("GSE211704/GSE211704_pheatmap_no_cluster.pdf",width=5,height=5)
print(p_list_GSE211704[[2]])
dev.off()
p_list_GSE127857=pheatmap_clincial(GSE127857_score,network1,GSE127857)
pdf("GSE127857/GSE127857_pheatmap.pdf",width=5,height=5)
print(p_list_GSE127857[[1]])
dev.off()
pdf("GSE127857/GSE127857_pheatmap_no_cluster.pdf",width=5,height=5)
print(p_list_GSE127857[[2]])
dev.off()
p_list_GSE207846=pheatmap_clincial(GSE207846_score,network1,GSE207846)
pdf("GSE207846/GSE207846_pheatmap.pdf",width=5,height=5)
print(p_list_GSE207846[[1]])
dev.off()
pdf("GSE207846/GSE207846_pheatmap_no_cluster.pdf",width=5,height=5)
print(p_list_GSE207846[[2]])
dev.off()

p_list_GSE164988=pheatmap_clincial(GSE164988_score,network1,GSE164988)
pdf("GSE164988/GSE164988_pheatmap.pdf",width=5,height=5)
print(p_list_GSE164988[[1]])
dev.off()
pdf("GSE164988/GSE1649886_pheatmap_no_cluster.pdf",width=5,height=5)
print(p_list_GSE164988[[2]])
dev.off()

p_list_combat=pheatmap_clincial(combat_score,network1,combat)
pdf("combat/combat_pheatmap.pdf",width=5,height=5)
print(p_list_combat[[1]])
dev.off()
pdf("combat/combat_pheatmap_no_cluster.pdf",width=5,height=5)
print(p_list_combat[[2]])
dev.off()
#################
##----重要的转录因子热图------
rm(list=ls())
library(magrittr)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
TF_path="D:/project/07.STAD_subtype/analysis0507_v7/04.OS/TF"
TF=read.table(file.path(TF_path,"single_cox_sig.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
TF_sig=unique(TF$Characteristics)
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA"
data=read.table(file.path(data_path,"TCGA-STAD-mRNA-tpm_log_case.txt"),header=TRUE,row.names=1, quote="",stringsAsFactors=F,check.names=F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
data=data[TF$Characteristics,]
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
train=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
train$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", train$sample)
test$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", test$sample)
##--function---#
pheatmap_data <- function(data,cluster){
  sample=intersect(colnames(data),cluster$sample)
  cluster=cluster[cluster$sample %in% sample,]
  cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
  cluster=cluster[order(cluster$cluster),]
  data=data[,cluster$sample]
  scaled_expr=t(scale(t(data),
                      center = T, scale=T)) 
  ha1=HeatmapAnnotation(df=data.frame(
    cluster=cluster$cluster),
    col=list(cluster=c("IM"= "#61b3de", "EP"="#ffee93")))
  p1=Heatmap(scaled_expr,name="TF",top_annotation=ha1,
             cluster_rows=TRUE,cluster_columns=FALSE,show_row_names=TRUE,
             show_column_names=FALSE,
             row_names_gp = gpar(fontsize=8))
  p2=Heatmap(scaled_expr,name="TF",top_annotation=ha1,
             cluster_rows=TRUE,cluster_columns=TRUE,show_row_names=TRUE,
             show_column_names=FALSE,
             row_names_gp = gpar(fontsize=8))
  clustercolor=c("#61b3de","#ffee93") 
  names(clustercolor)=c("IM","EP") #类型颜色
  ann_colors=list(cluster=clustercolor)
  annotation_col=data.frame(cluster=cluster$cluster)
  row.names(annotation_col)=cluster$sample
  p3=pheatmap(as.matrix(data),cluster_rows = T,cluster_cols = T,
              color=colorRampPalette(c("navy","white","firebrick3"))(100),
              fontsize_number =6,
              show_colnames = F,border_color = NA,scale = "row",show_rownames =T,
              annotation_col = annotation_col,# annotation_row = annotation_row,
              annotation_colors = ann_colors)
  p4=pheatmap(data,cluster_rows =T,cluster_cols = F,
              color=colorRampPalette(c("navy","white","firebrick3"))(100),
              show_colnames = F,border_color = NA,scale = "row",show_rownames =T,
              fontsize_number =6,
              annotation_col = annotation_col,# annotation_row = annotation_row,
              annotation_colors = ann_colors)
  p_list <- list(p1,p2,p3,p4)
  return(p_list)
}
setwd("D:/project/07.STAD_subtype/analysis0507_v7/18.heatmap")
dir.create("TF")
dir.create("TF/train")
dir.create("TF/test")
train_TF <- pheatmap_data(data,train)
pdf("TF/train/complexheatmap_no_cluster.pdf",width = 5,height=7.5)
print(train_TF[[1]])
dev.off()
pdf("TF/train/complexheatmap_cluster.pdf",width = 5,height=7.5)
print(train_TF[[2]])
dev.off()
pdf("TF/train/pheatmap_cluster.pdf",width = 5,height=7.5)
print(train_TF[[3]])
dev.off()
pdf("TF/train/pheatmap_no_cluster.pdf",width = 5,height=7.5)
print(train_TF[[4]])
dev.off()
test_TF <- pheatmap_data(data,test)
pdf("TF/test/complexheatmap_no_cluster.pdf",width = 5,height=7.5)
print(test_TF[[1]])
dev.off()
pdf("TF/test/complexheatmap_cluster.pdf",width = 5,height=7.5)
print(test_TF[[2]])
dev.off()
pdf("TF/test/pheatmap_cluster.pdf",width = 5,height=7.5)
print(test_TF[[3]])
dev.off()
pdf("TF/test/pheatmap_no_cluster.pdf",width = 5,height=7.5)
print(test_TF[[4]])
dev.off()
################
##----重要的甲基化位点热图------
rm(list=ls())
library(magrittr)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
meth_path="D:/project/07.STAD_subtype/analysis0507_v7/04.OS/CpG"
meth=read.table(file.path(meth_path,"single_cox_sig.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
meth_sig=unique(meth$Characteristics)
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/meth"
data=readRDS(file.path(data_path,"TCGA_GSE40279_promoter_knn_combat.rds"))
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
data=data[meth$Characteristics,]
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
train=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
train$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", train$sample)
test$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", test$sample)
##--function---#
pheatmap_data <- function(data,cluster){
  sample=intersect(colnames(data),cluster$sample)
  cluster=cluster[cluster$sample %in% sample,]
  cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
  cluster=cluster[order(cluster$cluster),]
  data=data[,cluster$sample]
  scaled_expr=t(scale(t(data),
                      center = T, scale=T)) 
  ha1=HeatmapAnnotation(df=data.frame(
    cluster=cluster$cluster),
    col=list(cluster=c("IM"= "#61b3de", "EP"="#ffee93")))
  p1=Heatmap(scaled_expr,name="CpG",top_annotation=ha1,
             cluster_rows=TRUE,cluster_columns=FALSE,show_row_names=FALSE,
             show_column_names=FALSE,
             row_names_gp = gpar(fontsize=8))
  p2=Heatmap(scaled_expr,name="CpG",top_annotation=ha1,
             cluster_rows=TRUE,cluster_columns=TRUE,show_row_names=FALSE,
             show_column_names=FALSE,
             row_names_gp = gpar(fontsize=8))
  clustercolor=c("#61b3de","#ffee93") 
  names(clustercolor)=c("IM","EP") #类型颜色
  ann_colors=list(cluster=clustercolor)
  annotation_col=data.frame(cluster=cluster$cluster)
  row.names(annotation_col)=cluster$sample
  p3=pheatmap(as.matrix(data),cluster_rows = T,cluster_cols = T,
              color=colorRampPalette(c("navy","white","firebrick3"))(100),
              fontsize_number =6,
              show_colnames = F,border_color = NA,scale = "row",show_rownames =F,
              annotation_col = annotation_col,# annotation_row = annotation_row,
              annotation_colors = ann_colors)
  p4=pheatmap(data,cluster_rows =T,cluster_cols = F,
              color=colorRampPalette(c("navy","white","firebrick3"))(100),
              show_colnames = F,border_color = NA,scale = "row",show_rownames =F,
              fontsize_number =6,
              annotation_col = annotation_col,# annotation_row = annotation_row,
              annotation_colors = ann_colors)
  p_list <- list(p1,p2,p3,p4)
  return(p_list)
}
setwd("D:/project/07.STAD_subtype/analysis0507_v7/18.heatmap")
dir.create("meth")
dir.create("meth/train")
dir.create("meth/test")
train_meth <- pheatmap_data(data,train)
pdf("meth/train/complexheatmap_no_cluster.pdf",width = 5,height=7.5)
print(train_meth[[1]])
dev.off()
pdf("meth/train/complexheatmap_cluster.pdf",width = 5,height=7.5)
print(train_meth[[2]])
dev.off()
pdf("meth/train/pheatmap_cluster.pdf",width = 5,height=7.5)
print(train_meth[[3]])
dev.off()
pdf("meth/train/pheatmap_no_cluster.pdf",width = 5,height=7.5)
print(train_meth[[4]])
dev.off()
test_meth <- pheatmap_data(data,test)
pdf("meth/test/complexheatmap_no_cluster.pdf",width = 5,height=7.5)
print(test_meth[[1]])
dev.off()
pdf("meth/test/complexheatmap_cluster.pdf",width = 5,height=7.5)
print(test_meth[[2]])
dev.off()
pdf("meth/test/pheatmap_cluster.pdf",width = 5,height=7.5)
print(test_meth[[3]])
dev.off()
pdf("meth/test/pheatmap_no_cluster.pdf",width = 5,height=7.5)
print(test_meth[[4]])
dev.off()
########
#################
##----重要的甲基化位点对应的基因的热图------
rm(list=ls())
library(magrittr)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
meth_path="D:/project/07.STAD_subtype/analysis0507_v7/02.meth_diff"
CpG_path="D:/project/07.STAD_subtype/analysis0507_v7/04.OS/CpG"
CpG=read.table(file.path(CpG_path,"single_cox_sig.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
CpG_gene=read.table(file.path(meth_path,"sig_cg_gene_point.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
CpG_gene=CpG_gene[CpG_gene$cg %in% CpG$Characteristics,]
gene=unique(CpG_gene$gene)
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA"
data=read.table(file.path(data_path,"TCGA-STAD-mRNA-tpm_log_case.txt"),header=TRUE,row.names=1, quote="",stringsAsFactors=F,check.names=F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
data=data[gene,]
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
train=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
train$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", train$sample)
test$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", test$sample)
##--function---#
pheatmap_data <- function(data,cluster){
  sample=intersect(colnames(data),cluster$sample)
  cluster=cluster[cluster$sample %in% sample,]
  cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
  cluster=cluster[order(cluster$cluster),]
  data=data[,cluster$sample]
  scaled_expr=t(scale(t(data),
                      center = T, scale=T)) 
  ha1=HeatmapAnnotation(df=data.frame(
    cluster=cluster$cluster),
    col=list(cluster=c("IM"= "#61b3de", "EP"="#ffee93")))
  p1=Heatmap(scaled_expr,name="CpG_gene",top_annotation=ha1,
             cluster_rows=TRUE,cluster_columns=FALSE,show_row_names=TRUE,
             show_column_names=FALSE,
             row_names_gp = gpar(fontsize=8))
  p2=Heatmap(scaled_expr,name="CpG_gene",top_annotation=ha1,
             cluster_rows=TRUE,cluster_columns=TRUE,show_row_names=TRUE,
             show_column_names=FALSE,
             row_names_gp = gpar(fontsize=8))
  clustercolor=c("#61b3de","#ffee93") 
  names(clustercolor)=c("IM","EP") #类型颜色
  ann_colors=list(cluster=clustercolor)
  annotation_col=data.frame(cluster=cluster$cluster)
  row.names(annotation_col)=cluster$sample
  p3=pheatmap(as.matrix(data),cluster_rows = T,cluster_cols = T,
              color=colorRampPalette(c("navy","white","firebrick3"))(100),
              fontsize_number =6,
              show_colnames = F,border_color = NA,scale = "row",show_rownames =T,
              annotation_col = annotation_col,# annotation_row = annotation_row,
              annotation_colors = ann_colors)
  p4=pheatmap(data,cluster_rows =T,cluster_cols = F,
              color=colorRampPalette(c("navy","white","firebrick3"))(100),
              show_colnames = F,border_color = NA,scale = "row",show_rownames =T,
              fontsize_number =6,
              annotation_col = annotation_col,# annotation_row = annotation_row,
              annotation_colors = ann_colors)
  p_list <- list(p1,p2,p3,p4)
  return(p_list)
}
setwd("D:/project/07.STAD_subtype/analysis0507_v7/18.heatmap")
dir.create("CpG_gene")
dir.create("CpG_gene/train")
dir.create("CpG_gene/test")
train_gene <- pheatmap_data(data,train)
pdf("CpG_gene/train/complexheatmap_no_cluster.pdf",width = 5,height=7.5)
print(train_gene[[1]])
dev.off()
pdf("CpG_gene/train/complexheatmap_cluster.pdf",width = 5,height=7.5)
print(train_gene[[2]])
dev.off()
pdf("CpG_gene/train/pheatmap_cluster.pdf",width = 5,height=7.5)
print(train_gene[[3]])
dev.off()
pdf("CpG_gene/train/pheatmap_no_cluster.pdf",width = 5,height=7.5)
print(train_gene[[4]])
dev.off()
test_gene <- pheatmap_data(data,test)
pdf("CpG_gene/test/complexheatmap_no_cluster.pdf",width = 5,height=7.5)
print(test_gene[[1]])
dev.off()
pdf("CpG_gene/test/complexheatmap_cluster.pdf",width = 5,height=7.5)
print(test_gene[[2]])
dev.off()
pdf("CpG_gene/test/pheatmap_cluster.pdf",width = 5,height=7.5)
print(test_gene[[3]])
dev.off()
pdf("CpG_gene/test/pheatmap_no_cluster.pdf",width = 5,height=7.5)
print(test_gene[[4]])
dev.off()
################
##--GEO数据的甲基化位点情况---#
###--GEO 重要的甲基化位点热图-------
rm(list=ls())
library(magrittr)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
meth_path="D:/project/07.STAD_subtype/analysis0507_v7/04.OS/CpG"
meth=read.table(file.path(meth_path,"single_cox_sig.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
meth_sig=unique(meth$Characteristics)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/17.cluster_valiation")
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/GEO_methylation"
GSE211704_data=read.table(file.path(data_path,"GSE211704/GSE211704_Tumor_knn.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
GSE211704_data=GSE211704_data[meth_sig,]
GSE211704=read.table("GSE211704/ConsensusClusterPlus10-100-km-euclidean/consensusClass2.txt",sep="\t",header=F,stringsAsFactors = F,quote = "")
colnames(GSE211704)=c("sample","cluster")
GSE211704$cluster=paste0("cluster",GSE211704$cluster)
GSE211704$cluster[GSE211704$cluster=="cluster1"]="IM"
GSE211704$cluster[GSE211704$cluster=="cluster2"]="EP"
GSE127857_data=read.table(file.path(data_path,"GSE127857/GSE127857_promoter.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
GSE127857_data=GSE127857_data[meth_sig,]
GSE127857=read.table("GSE127857/ConsensusClusterPlus10-100-km-euclidean/consensusClass2.txt",sep="\t",header=F,stringsAsFactors = F,quote = "")
colnames(GSE127857)=c("sample","cluster")
GSE127857$cluster=paste0("cluster",GSE127857$cluster)
GSE127857$cluster[GSE127857$cluster=="cluster2"]="IM"
GSE127857$cluster[GSE127857$cluster=="cluster1"]="EP"
GSE207846_data=read.table(file.path(data_path,"GSE207846/GSE207846_Tumor_knn.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
GSE207846_data=GSE207846_data[meth_sig,]
GSE207846=read.table("GSE207846/ConsensusClusterPlus10-100-km-euclidean/consensusClass2.txt",sep="\t",header=F,stringsAsFactors = F,quote = "")
colnames(GSE207846)=c("sample","cluster")
GSE207846$cluster=paste0("cluster",GSE207846$cluster)
GSE207846$cluster[GSE207846$cluster=="cluster2"]="IM"
GSE207846$cluster[GSE207846$cluster=="cluster1"]="EP"
GSE164988_data=read.table(file.path(data_path,"GSE164988/GSE164988_knn.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
GSE164988_data=GSE164988_data[meth_sig,]
GSE164988=read.table("GSE164988/ConsensusClusterPlus10-100-km-euclidean/consensusClass2.txt",sep="\t",header=F,stringsAsFactors = F,quote = "")
colnames(GSE164988)=c("sample","cluster")
GSE164988$cluster=paste0("cluster",GSE164988$cluster)
GSE164988$cluster[GSE164988$cluster=="cluster1"]="IM"
GSE164988$cluster[GSE164988$cluster=="cluster2"]="EP"
combat=read.table("combat/ConsensusClusterPlus10-100-km-euclidean/consensusClass2.txt",sep="\t",header=F,stringsAsFactors = F,quote = "")
colnames(combat)=c("sample","cluster")
combat$cluster=paste0("cluster",combat$cluster)
combat$cluster[combat$cluster=="cluster1"]="IM"
combat$cluster[combat$cluster=="cluster2"]="EP"
combat_data=read.table(file.path(data_path,"combat/combat_data.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
combat_data=combat_data[meth_sig,]
pheatmap_data <- function(data,cluster){
  sample=intersect(colnames(data),cluster$sample)
  cluster=cluster[cluster$sample %in% sample,]
  cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
  cluster=cluster[order(cluster$cluster),]
  data=data[,cluster$sample]
  scaled_expr=t(scale(t(data),
                      center = T, scale=T)) 
  ha1=HeatmapAnnotation(df=data.frame(
    cluster=cluster$cluster),
    col=list(cluster=c("IM"= "#61b3de", "EP"="#ffee93")))
  p1=Heatmap(scaled_expr,name="CpG",top_annotation=ha1,
             cluster_rows=TRUE,cluster_columns=FALSE,show_row_names=FALSE,
             show_column_names=FALSE,
             row_names_gp = gpar(fontsize=8))
  p2=Heatmap(scaled_expr,name="CpG",top_annotation=ha1,
             cluster_rows=TRUE,cluster_columns=TRUE,show_row_names=FALSE,
             show_column_names=FALSE,
             row_names_gp = gpar(fontsize=8))
  clustercolor=c("#61b3de","#ffee93") 
  names(clustercolor)=c("IM","EP") #类型颜色
  ann_colors=list(cluster=clustercolor)
  annotation_col=data.frame(cluster=cluster$cluster)
  row.names(annotation_col)=cluster$sample
  p3=pheatmap(as.matrix(data),cluster_rows = T,cluster_cols = T,
              color=colorRampPalette(c("navy","white","firebrick3"))(100),
              fontsize_number =6,
              show_colnames = F,border_color = NA,scale = "row",show_rownames =F,
              annotation_col = annotation_col,# annotation_row = annotation_row,
              annotation_colors = ann_colors)
  p4=pheatmap(data,cluster_rows =T,cluster_cols = F,
              color=colorRampPalette(c("navy","white","firebrick3"))(100),
              show_colnames = F,border_color = NA,scale = "row",show_rownames =F,
              fontsize_number =6,
              annotation_col = annotation_col,# annotation_row = annotation_row,
              annotation_colors = ann_colors)
  p_list <- list(p1,p2,p3,p4)
  return(p_list)
}
setwd("D:/project/07.STAD_subtype/analysis0507_v7/18.heatmap/meth")
dir.create("GSE211704")
dir.create("GSE127857")
dir.create("GSE207846")
dir.create("GSE164988")
dir.create("combat")
p_list_GSE211704=pheatmap_data(GSE211704_data,GSE211704)
pdf("GSE211704/GSE211704_complexheatmap_no_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE211704[[1]])
dev.off()
pdf("GSE211704/GSE211704_complexheatmap_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE211704[[2]])
dev.off()
pdf("GSE211704/GSE211704_pheatmap_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE211704[[3]])
dev.off()
pdf("GSE211704/GSE211704_pheatmap_no_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE211704[[4]])
dev.off()
p_list_GSE127857=pheatmap_data(GSE127857_data,GSE127857)
pdf("GSE127857/GSE127857_complexheatmap_no_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE127857[[1]])
dev.off()
pdf("GSE127857/GSE127857_complexheatmap_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE127857[[2]])
dev.off()
pdf("GSE127857/GSE127857_pheatmap_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE127857[[3]])
dev.off()
pdf("GSE127857/GSE127857_pheatmap_no_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE127857[[4]])
dev.off()
p_list_GSE207846=pheatmap_data(GSE207846_data,GSE207846)
pdf("GSE207846/GSE207846_complexheatmap_no_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE207846[[1]])
dev.off()
pdf("GSE207846/GSE207846_complexheatmap_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE207846[[2]])
dev.off()
pdf("GSE207846/GSE207846_pheatmap_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE207846[[3]])
dev.off()
pdf("GSE207846/GSE207846_pheatmap_no_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE207846[[4]])
dev.off()
p_list_GSE164988=pheatmap_data(GSE164988_data,GSE164988)
pdf("GSE164988/GSE164988_complexheatmap_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE164988[[1]])
dev.off()
pdf("GSE164988/GSE1649886_complexheatmap_no_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE164988[[2]])
dev.off()
pdf("GSE164988/GSE164988_pheatmap_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE164988[[3]])
dev.off()
pdf("GSE164988/GSE164988_pheatmap_no_cluster.pdf",width = 5,height=7.5)
print(p_list_GSE164988[[4]])
dev.off()
p_list_combat=pheatmap_data(combat_data,combat)
pdf("combat/combat_complexheatmap_cluster.pdf",width = 5,height=7.5)
print(p_list_combat[[2]])
dev.off()
pdf("combat/combat_complexheatmap_no_cluster.pdf",width = 5,height=7.5)
print(p_list_combat[[1]])
dev.off()
pdf("combat/combat_pheatmap_cluster.pdf",width = 5,height=7.5)
print(p_list_combat[[3]])
dev.off()
pdf("combat/combat_pheatmap_no_cluster.pdf",width = 5,height=7.5)
print(p_list_combat[[4]])
dev.off()
###