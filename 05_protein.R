#################
##--下载蛋白质数据----
rm(list=ls())
library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
library(reshape2)
query <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification"
)
GDCdownload(query)
GDCprepare(query, save = T,save.filename = "TCGA-STAD_protein.Rdata")
#########
#----使用k近邻法对蛋白质NA的补缺失值------
rm(list=ls())
library(DMwR2)
protein_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/protein"
load(file.path(protein_path,"TCGA-STAD_protein.Rdata"))
data=as.data.frame(data)
rownames(data)=data$peptide_target
data=data[,-c(1:5)]
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
samplesNT=colnames(data)[which(group==1)]
samplesTP=colnames(data)[which(group==0)]
type=data.frame(sample=c(samplesNT,samplesTP),type=c(rep("Normal",length(samplesNT)),
                                                     rep("Tumor",length(samplesTP))))
del_row=which(rowSums(is.na(data))/ncol(data) > 0.5)
data=data[-del_row, ]
data1=knnImputation(data)
anyNA(data1)
saveRDS(data1,file.path(protein_path,"TCGA-STAD_protein_knn.rds"))
##都是tumor
#######
##---蛋白质差异分析----------
rm(list=ls())
library(reshape2)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/")
dir.create("13.protein")
setwd("13.protein")
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
cluster=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", cluster$sample)
rownames(cluster)=cluster$sample
cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
protein_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/protein"
data=readRDS(file.path(protein_path,"TCGA-STAD_protein_knn.rds"))
sample=intersect(colnames(data),cluster$sample)
cluster=cluster[cluster$sample %in% sample,]
data_case=data[,cluster$sample[cluster$cluster=="IM"]]
data_control=data[,cluster$sample[cluster$cluster=="EP"]]
mean_case=apply(data_case,1,mean)
mean_control=apply(data_control,1,mean)
de=mean_case-mean_control
fold_change=abs(mean_case/mean_control)
loc=which(fold_change>=2|fold_change<=1/2)
data=cbind(data_case,data_control)
data$de=de
data$fc=fold_change
#############t检验
p=NULL
sig_fold_change=which(fold_change>=2|fold_change<=1/2)#
for(i in 1:dim(data)[1]){
  p_value=t.test(data_case[i,],data_control[i,])$p.value
  p=c(p,p_value)
}
data$p=p
data$padj=p.adjust(data$p,method="fdr",length(data$p))
sig_padj=which(data$padj<0.05)#
sig_m=intersect(sig_fold_change,sig_padj)###
fc.padj=data.frame(id=rownames(data),fc=data$fc,p=data$p,padj=data$padj)
sig.data=data[sig_m,]#
sig.protein=rownames(sig.data)
protein_result=list(data,fc.padj,sig.data,sig.protein)
saveRDS(protein_result,"protein_diff_result.rds")
write.table(sig.data,"sig_protein.xls",sep="\t",quote=FALSE,row.names = T)

##---画热图-------
rm(list=ls())
library(magrittr)
library(ComplexHeatmap)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/13.protein")
data=read.delim("sig_protein.xls",sep="\t",header=T,stringsAsFactors=F,check.names=F,quote="")
data$fc=abs(data$fc)
data=data[order(data$fc),]
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
cluster=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", cluster$sample)
rownames(cluster)=cluster$sample
cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
sample=intersect(colnames(data),cluster$sample)
cluster=cluster[cluster$sample %in% sample,]
cluster=cluster[order(cluster$cluster),]
Up=data%>%dplyr::filter(de>0)
Down=data%>%dplyr::filter(de<0)
expr=rbind(Up,Down)
expr=expr[,cluster$sample]
###---标准化----###
###---标准化----###
scaled_expr=t(scale(t(expr))) 
ha1=HeatmapAnnotation(df=data.frame(
  cluster=cluster$cluster),
  col=list(cluster=c("IM"= "#61b3de", "EP"="#ffee93")))
protein=c(rep("Up",nrow(Up)),rep("Down",nrow(Down)))
ha2=rowAnnotation(df=data.frame(
  protein=c(rep("Up",nrow(Up)),rep("Down",nrow(Down)))),
  col=list(protein=c("Up"="tomato","Down"="steelblue")))
library(circlize)
pdf(file="complexheatmap_nocluster.pdf",width=7,height=7)
Heatmap(scaled_expr,name="expression",top_annotation=ha1,left_annotation =ha2 ,
        cluster_rows=TRUE,cluster_columns=FALSE,show_row_names=TRUE,
        show_column_names=FALSE,
        column_split = cluster$cluster,
        row_split = protein,
        row_names_gp = gpar(fontsize=8),
        column_names_gp=gpar(fontsize=8))
dev.off()

library(pheatmap)
annotation_col = data.frame( cluster = cluster$cluster)
rownames(annotation_col) = cluster$sample
annotation_row = data.frame( protein=c(rep("Up",nrow(Up)),rep("Down",nrow(Down))))
rownames(annotation_row) = c(rownames(Up),rownames(Down))
ann_colors = list( cluster = c(IM ="#61b3de", EP = "#ffee93"),
                   protein=c(Up="#61b3de",Down="#ffee93"))
# scaled_expr[scaled_expr>5]=5
# scaled_expr[scaled_expr< -5]= -5
bk <- c(seq(-5,5,by=0.1)) #定义热图中值的区间
cluster1.length=length(which(cluster$cluster=="IM"))
up_length=dim(Up)[1]
pdf(file="pheatmap.pdf",width=5,height=6.5)
pheatmap(expr,
         scale = "row",
         cluster_row = F, cluster_col = F, border=NA,
         #display_numbers = pmt,
         show_colnames = F,
         show_rownames = T,
         color = colorRampPalette(rev(c("#BC3C29FF","#FCFCE3","#0072B5FF")))(length(bk)),
         fontsize_number =8, number_color = "white",
         gaps_col =cluster1.length,
         gaps_row = up_length,
         annotation_row = annotation_row,
         
         annotation_col =annotation_col, 
         annotation_colors =ann_colors)
dev.off()