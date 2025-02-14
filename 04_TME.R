##---免疫细胞比较 两类-------
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(ggprism)
library(ggsci)
library(cowplot)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/07.scRNA-seq/harmony_0.5_30/markers_by_self")
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
train=read.table(file.path(cluster_path,"train_cluster.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table(file.path(cluster_path,"test_cluster.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
data=read.table('ssGSEA/gsva_matrix.txt',sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
test$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", test$sample)
train$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", train$sample)
data=cbind(cell=rownames(data),data)
data1=data %>%
  gather(sample, score, -cell)
ssGSEA_comapre <- function(data1,cluster,title1,cells){
  data=merge(data1,cluster,by="sample")
  plot.data=data
  plot.data$cluster=factor(plot.data$cluster,c("cluster1","cluster2"))
  plot.data$cell=factor(plot.data$cell,levels=cells)
  p=ggboxplot(plot.data, x="cell", y="score", color="cluster",
              xlab="",
              ylab="score",
              title=paste0("scRNA-seq Deconvolution of ",title1),
              legend.title="",
              width=0.8,
              palette = c("#61b3de", "#ffee93"))+
    rotate_x_text(90)+ ##这个是改倾斜度的
    stat_compare_means(aes(group=cluster),
                       method="wilcox.test",
                       symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
  return(print(p))
}
cells=unique(data1$cell)
p_train=ssGSEA_comapre(data1,train,"Train Set",cells)
p_test=ssGSEA_comapre(data1,test,"Test Set",cells)
ssGSEA_list=list(p_train,p_test)
pdf("ssGSEA/ssGSEA_cluster.pdf",width=10,height=5.5)
plot_grid(plotlist=ssGSEA_list,ncol=2)
dev.off()
#################
###------以下全是两类的分析-----------
##---根据单细胞的结果 给cluster重新定义名字----
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2")
train=read.table("train_cluster.xls",sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table("test_cluster.xls",sep="\t",header=T,stringsAsFactors = F,quote="")
train$cluster[train$cluster=="cluster1"]="IM"
train$cluster[train$cluster=="cluster2"]="EP"
test$cluster[test$cluster=="cluster1"]="IM"
test$cluster[test$cluster=="cluster2"]="EP"
write.table(train,"train_cluster_rename.xls",sep="\t",quote=FALSE,row.names = F)
write.table(test,"test_cluster_rename.xls",sep="\t",quote=FALSE,row.names = F)
##---survival 重命名----
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/")
train=read.table("ConsensusClusterPlus10-100-km-euclidean/survival.txt",header=T,sep="\t",stringsAsFactors = F)
train$class[train$class==2]="IM"
train$class[train$class==1]="EP"
title_train="Train Set"
##---test survival---##
test=read.table("test_ConsensusClusterPlus10-100-km-euclidean/test_survival.txt",header=T,sep="\t",stringsAsFactors = F)
test$class[test$class==1]="IM"
test$class[test$class==2]="EP"
write.table(train,"cluster_divide2/train_survival_rename.xls",sep="\t",quote=FALSE,row.names = F)
write.table(test,"cluster_divide2/test_survival_rename.xls",sep="\t",quote=FALSE,row.names = F)
##------十、免疫微环境比较------------
##----多种免疫细胞浸润分析----------
#install.packages("devtools")
#install_github("icbi-lab/immunedeconv")
#https://icbi-lab.github.io/immunedeconv/articles/immunedeconv.html
#install.packages("remotes")
#remotes::install_github("icbi-lab/immunedeconv")
#####癌症样本
#####xcell-----
rm(list=ls())
library(immunedeconv)
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA"
exprMatrix = read.table(file.path(data_path,"TCGA-STAD-mRNA-tpm_log_case.txt"),header=TRUE,row.names=1, quote="",stringsAsFactors=F,check.names=F)
res_xcell = deconvolute(exprMatrix, method="xcell")
setwd("D:/project/07.STAD_subtype/analysis0507_v7/08.immune_cell")
dir.create("xcell")
write.table(res_xcell, file="xcell/xcell.txt", sep="\t", col.names=T, row.names=F, quote=F)
#####epic------
res_epic = deconvolute(exprMatrix,method= "epic")
dir.create("epic")
write.table(res_epic, file="epic/epic.txt", sep="\t", col.names=T, row.names=F, quote=F)
#####quantiseq-----
res_quantiseq = deconvolute(exprMatrix, method="quantiseq")
dir.create("quantiseq")
write.table(res_quantiseq, file="quantiseq/quantiseq.txt", sep="\t", col.names=T, row.names=F, quote=F)
#####timer----
res_timer = deconvolute(exprMatrix, method="timer",indications=rep("STAD",dim(exprMatrix)[2])) #每个样本的癌症类型都要指出
dir.create("timer")
write.table(res_timer, file="timer/timer.txt", sep="\t", col.names=T, row.names=F, quote=F)
#####CIBERSORT----
source('D:/project/07.STAD_subtype/analysis0507_v7/08.immune_cell/Cibersort.R')
res_cibersort=CIBERSORT('D:/project/07.STAD_subtype/analysis0507_v7/08.immune_cell/LM22.txt',
                        'D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA/TCGA-STAD-mRNA-tpm_log_case.txt', perm = 100, QN = F) #perm置换次数=1000，QN分位数归一化
dir.create("cibersort")
write.table(res_cibersort, file="cibersort/cibersort.txt", sep="\t", col.names=T, row.names=T, quote=F)
####ssGSEA-------
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA"
expression = read.table(file.path(data_path,"TCGA-STAD-mRNA-tpm_log_case.txt"),header=TRUE,row.names=1, quote="",stringsAsFactors=F,check.names=F)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/08.immune_cell")
##读取背景基因集合
gene_set=read.table("gene_set.txt",header=T,sep="\t",stringsAsFactors = F)[,c(1:2)]
list=split(as.matrix(gene_set)[,1], gene_set[,2])##存储的是每个免疫细胞对应的基因,构建背景基因集合
library(GSVA)
gsva_matrix=gsva(as.matrix(expression),list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)##进行gsva分析
dir.create("ssGSEA")
write.table(gsva_matrix,"ssGSEA/gsva_matrix.txt",sep="\t",quote=FALSE,row.names = TRUE)##
####ssGSEA immport----
rm(list=ls())
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA"
expression = read.table(file.path(data_path,"TCGA-STAD-mRNA-tpm_log_case.txt"),header=TRUE,row.names=1, quote="",stringsAsFactors=F,check.names=F)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/08.immune_cell")
##读取背景基因集合
gene_set=read.table("immune_reaction_genesets.txt",header=T,sep="\t",stringsAsFactors = F)[,c(1:2)]
list=split(as.matrix(gene_set)[,1], gene_set[,2])##存储的是每个免疫细胞对应的基因,构建背景基因集合
library(GSVA)
gsva_matrix=gsva(as.matrix(expression),list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)##进行gsva分析
dir.create("immport")
write.table(gsva_matrix,"immport/immune_reaction_gsva_matrix.txt",sep="\t",quote=FALSE,row.names = TRUE)##输出结果
###################
##--estimate-----
rm(list=ls())
library(limma)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/08.immune_cell")
dir.create("estimate")
setwd("estimate")
#运行estimate包
filterCommonGenes(input.f="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA/TCGA-STAD-mRNA-tpm_log_case.txt", ##表达数据
                  output.f="ESTIMATE_input.gct", ## #生成ESTIMATE 的输入文件
                  id="GeneSymbol")
# 该功能将每个平台的不同数量的基因与10412个普通基因相结合
#这个功能计算基质，免疫，并估计得分每个样本使用基因表达数据
estimateScore("ESTIMATE_input.gct", "ESTIMATE_score.gct")
plotPurity(scores="ESTIMATE_score.gct", samples="all_samples")
#根据ESTIMATE score绘制肿瘤纯度
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/08.immune_cell/estimate")
#将评分保存为txt格式
ESTIMATE_score=read.table("ESTIMATE_score.gct", skip = 2,header = TRUE,row.names = 1,check.names=F)#前两行跳过
colnames(ESTIMATE_score)=gsub("[.]","-",colnames(ESTIMATE_score))
ESTIMATE_score=ESTIMATE_score[,2:ncol(ESTIMATE_score)]
write.table(ESTIMATE_score,file = "ESTIMATE_score.txt",quote = F,sep = "\t",row.names = T)
##---estimate comapre ----
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(ggprism)
library(ggsci)
library(cowplot)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/08.immune_cell")
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
train=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
setwd("D:/project/07.STAD_subtype/analysis0507_v7/08.immune_cell/estimate")
data=read.table("ESTIMATE_score.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names=F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
test$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", test$sample)
train$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", train$sample)
data=cbind(type=rownames(data),data)
data=data[-4,]
data1=data %>%
  gather(sample, score, -type)
estimate_comapre <- function(data1,cluster,title1){
  data=merge(data1,cluster,by="sample")
  data1=data[,-1]
  colnames(data1)=c("type","score","class")
  data1$class=factor(data1$class,c("IM","EP"))
  p=ggviolin(data1, x="type", y="score", color="class",
             xlab="",
             ylab="Score",
             title=paste0("Estimate of ",title1),
             legend.title="",
             width=0.8,
             palette = c("#61b3de", "#ffee93"),add='boxplot')+
    rotate_x_text(30)+ ##这个是改倾斜度的
    stat_compare_means(aes(group=class),
                       method="wilcox.test",
                       symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
  return(print(p))
}
p_train=estimate_comapre(data1,train,"Train Set")
p_test=estimate_comapre(data1,test,"Test Set")
estimate_list=list(p_train,p_test)
pdf("estimate_cluster.pdf",width=8,height=4)
plot_grid(plotlist=estimate_list,ncol=2)
dev.off()
##---estimate compare half Violin----
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(ggprism)
library(ggsci)
library(cowplot)
library(ggplot2)
library(gghalves)
library(tidyverse)
# 绘图要调用的函数
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
setwd("D:/project/07.STAD_subtype/analysis0507_v7/08.immune_cell")
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
train=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
setwd("D:/project/07.STAD_subtype/analysis0507_v7/08.immune_cell/estimate")
data=read.table("ESTIMATE_score.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names=F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
test$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", test$sample)
train$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", train$sample)
data=cbind(type=rownames(data),data)
data=data[-4,]
data1=data %>%
  gather(sample, score, -type)
estimate_comapre <- function(data1,cluster,title1){
  data=merge(data1,cluster,by="sample")
  data1=data[,-1]
  colnames(data1)=c("type","score","class")
  data1$class=factor(data1$class,c("IM","EP"))
  p1=ggplot(data1, aes(x = type,y = score, fill = class))+
    geom_split_violin(trim = T,colour="white")+
    geom_point(stat = 'summary',fun=mean,
               position = position_dodge(width = 0.2))+
    scale_fill_manual(values = c("#61b3de", "#ffee93"))+
    stat_summary(fun.min = function(x){quantile(x)[2]},
                 fun.max = function(x){quantile(x)[4]},
                 geom = 'errorbar',color='black',
                 width=0.01,size=0.5,
                 position = position_dodge(width = 0.2))+
    stat_compare_means(data = data1, aes(x = type,y = score),
                       # 修改显著性标注：
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "-")),
                       label = "p.signif",
                       label.y = max(data1$score),
                       hide.ns = F)+
    theme_bw()+labs(y="Score",x="",title=paste0("Estimate of ",title1))+
    theme(axis.text.x = element_text(angle = 30, hjust = 1), 
          legend.position = "top",
          #legend.key = element_rect(fill = c("#1ba7b3","#dfb424")),
          legend.justification = "right")
  
}

p_train=estimate_comapre(data1,train,"Train Set")
p_test=estimate_comapre(data1,test,"Test Set")
estimate_list=list(p_train,p_test)
pdf("estimate_cluster_v2.pdf",width=6.5,height=5)
plot_grid(plotlist=estimate_list,ncol=2)
dev.off()
pdf("estimate_train_v2.pdf",width=3,height=5)
p_train
dev.off()
pdf("estimate_test_v2.pdf",width=3,height=5)
p_test
dev.off()




###########
##---免疫细胞比较 单独画 每个画箱式图--------
##---免疫细胞比较 train------
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(ggprism)
library(ggsci)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/08.immune_cell")
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
cluster=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", cluster$sample)
rownames(cluster)=cluster$sample
##--ssGSEA train------
data=read.table('ssGSEA/gsva_matrix.txt',sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
data=cbind(cell=rownames(data),data)
data1=data %>%
  gather(sample, score, -cell)
data=merge(data1,cluster,by="sample")
plot.data=data
plot.data$cluster=factor(plot.data$cluster,c("IM","EP"))
p=ggboxplot(plot.data, x="cell", y="score", color="cluster",
            xlab="",
            ylab="score",
            title="ssGSEA",
            legend.title="",
            width=0.8,
            palette = c("#61b3de", "#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=cluster),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("ssGSEA/compare_train.pdf",width=10,height=6)
p
dev.off()
##---epic train-----
data=read.table("epic/epic.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
data=data[-which(data$cell_type=="uncharacterized cell"),]
data1=data %>%
  gather(sample, score, -cell_type)
data=merge(data1,cluster,by="sample")
data1=data[,-1]
colnames(data1)=c("cell","fraction","class")
data1$class=factor(data1$class,c("IM","EP"))
p=ggboxplot(data1, x="cell", y="fraction", color="class",
            xlab="",
            ylab="Fraction",
            title="epic",
            legend.title="Type",
            width=0.8,
            palette = c("#61b3de", "#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=class),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("epic/cluster_train.pdf",width=8,height=6)
p
dev.off()
##--quantiseq train--------
data=read.table("quantiseq/quantiseq.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
data=data[-which(data$cell_type=="uncharacterized cell"),]
data1=data %>%
  gather(sample, score, -cell_type)
data=merge(data1,cluster,by="sample")
data1=data[,-1]
colnames(data1)=c("cell","fraction","class")
data1$class=factor(data1$class,c("IM","EP"))
p=ggboxplot(data1, x="cell", y="fraction", color="class",
            xlab="",
            ylab="Fraction",
            title="quantiseq",
            legend.title="Type",
            width=0.8,
            palette = c("#61b3de", "#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=class),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("quantiseq/cluster_train.pdf",width=8,height=6)
p
dev.off()
##----xcell train----------
data=read.table("xcell/xcell.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
data1=data %>%
  gather(sample, score, -cell_type)
data=merge(data1,cluster,by="sample")
data1=data[,-1]
colnames(data1)=c("cell","fraction","class")
data1$class=factor(data1$class,c("IM","EP"))
p=ggboxplot(data1, x="cell", y="fraction", color="class",
            xlab="",
            ylab="Fraction",
            title="xcell",
            legend.title="Type",
            width=0.8,
            palette = c("#61b3de", "#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=class),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("xcell/cluster_train.pdf",width=10,height=6)
p
dev.off()


##---timer train----
data=read.table("timer/timer.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
data1=data %>%
  gather(sample, score, -cell_type)
data=merge(data1,cluster,by="sample")
data1=data[,-1]
colnames(data1)=c("cell","fraction","class")
data1$class=factor(data1$class,c("IM","EP"))
p=ggboxplot(data1, x="cell", y="fraction", color="class",
            xlab="",
            ylab="Fraction",
            title="timer",
            legend.title="Type",
            width=0.8,
            palette = c("#61b3de", "#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=class),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("timer/cluster_train.pdf",width=8,height=6)
p
dev.off()


###---cibersort得到的是占比------
data=read.table("cibersort/cibersort.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
data=data[data[,"P-value"]<0.05,]
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", rownames(data))
cluster1=cluster
rownames(cluster1)=cluster1$sample
samesample=intersect(rownames(data),cluster1$sample)
cluster1=cluster1[samesample,]
data=data[samesample,]
result=NULL
for(i in 1:22){
  result.sub=data.frame(cell=rep(colnames(data)[i],dim(data)[1]),
                        fraction=data[,i],
                        type=cluster1$cluster)
  colnames(result.sub)=c("cell","fraction","class")
  result=rbind(result,result.sub)
}
library(ggpubr)
plot.data=result
library(ggplot2)
library(ggprism)
library(ggsci)
plot.data$class=factor(plot.data$class,c("IM","EP"))
p=ggboxplot(plot.data, x="cell", y="fraction", color="class",
            xlab="",
            ylab="Fraction",
            title="Cibersort",
            legend.title="Type",
            width=0.8,
            palette = c("#61b3de", "#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=class),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("Cibersort/compare_train.pdf",width=10,height=6)
p
dev.off()
##---immport----
data=read.table('immport/immune_reaction_gsva_matrix.txt',sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
data=cbind(pathway=rownames(data),data)
data1=data %>%
  gather(sample, score, -pathway)
data=merge(data1,cluster,by="sample")
plot.data=data
plot.data$cluster=factor(plot.data$cluster,c("IM","EP"))
p=ggboxplot(plot.data, x="pathway", y="score", color="cluster",
            xlab="",
            ylab="score",
            title="ssGSEA",
            legend.title="",
            width=0.8,
            palette = c("#61b3de", "#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=cluster),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("immport/train_compare.pdf",width=10,height=6)
p
dev.off()
##################
###############
##---免疫细胞比较 test------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/08.immune_cell")
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
cluster=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", cluster$sample)
library(tidyverse)
library(ggpubr)
library(ggprism)
library(ggsci)
##--ssGSEA test ------
data=read.table('ssGSEA/gsva_matrix.txt',sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
data=cbind(cell=rownames(data),data)
data1=data %>%
  gather(sample, score, -cell)
data=merge(data1,cluster,by="sample")
plot.data=data
plot.data$cluster=factor(plot.data$cluster,c("IM","EP"))
p=ggboxplot(plot.data, x="cell", y="score", color="cluster",
            xlab="",
            ylab="score",
            title="ssGSEA test",
            legend.title="",
            width=0.8,
            palette = c("#61b3de", "#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=cluster),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("ssGSEA/compare_test.pdf",width=10,height=6)
p
dev.off()
##---epic test -----
data=read.table("epic/epic.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
data=data[-which(data$cell_type=="uncharacterized cell"),]
data1=data %>%
  gather(sample, score, -cell_type)
data=merge(data1,cluster,by="sample")
data1=data[,-1]
colnames(data1)=c("cell","fraction","class")
data1$class=factor(data1$class,c("IM","EP"))
p=ggboxplot(data1, x="cell", y="fraction", color="class",
            xlab="",
            ylab="Fraction",
            title="epic test",
            legend.title="Type",
            width=0.8,
            palette = c("#61b3de", "#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=class),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("epic/cluster_test.pdf",width=8,height=6)
p
dev.off()
##--quantiseq test--------
data=read.table("quantiseq/quantiseq.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
data=data[-which(data$cell_type=="uncharacterized cell"),]
data1=data %>%
  gather(sample, score, -cell_type)
data=merge(data1,cluster,by="sample")
data1=data[,-1]
colnames(data1)=c("cell","fraction","class")
data1$class=factor(data1$class,c("IM","EP"))
p=ggboxplot(data1, x="cell", y="fraction", color="class",
            xlab="",
            ylab="Fraction",
            title="quantiseq test",
            legend.title="Type",
            width=0.8,
            palette = c("#61b3de", "#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=class),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("quantiseq/cluster_test.pdf",width=8,height=6)
p
dev.off()
##----xcell test----------
data=read.table("xcell/xcell.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
data1=data %>%
  gather(sample, score, -cell_type)
data=merge(data1,cluster,by="sample")
data1=data[,-1]
colnames(data1)=c("cell","fraction","class")
data1$class=factor(data1$class,c("IM","EP"))
p=ggboxplot(data1, x="cell", y="fraction", color="class",
            xlab="",
            ylab="Fraction",
            title="xcell test",
            legend.title="Type",
            width=0.8,
            palette = c("#61b3de", "#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=class),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("xcell/cluster_test.pdf",width=10,height=6)
p
dev.off()


##---timer test----
data=read.table("timer/timer.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
data1=data %>%
  gather(sample, score, -cell_type)
data=merge(data1,cluster,by="sample")
data1=data[,-1]
colnames(data1)=c("cell","fraction","class")
data1$class=factor(data1$class,c("IM","EP"))
p=ggboxplot(data1, x="cell", y="fraction", color="class",
            xlab="",
            ylab="Fraction",
            title="timer test",
            legend.title="Type",
            width=0.8,
            palette = c("#61b3de", "#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=class),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("timer/cluster_test.pdf",width=8,height=6)
p
dev.off()


###---cibersort test得到的是占比------
data=read.table("cibersort/cibersort.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
data=data[data[,"P-value"]<0.05,]
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", rownames(data))
cluster1=cluster
rownames(cluster1)=cluster1$sample
samesample=intersect(rownames(data),cluster1$sample)
cluster1=cluster1[samesample,]
data=data[samesample,]
result=NULL
for(i in 1:22){
  result.sub=data.frame(cell=rep(colnames(data)[i],dim(data)[1]),
                        fraction=data[,i],
                        type=cluster1$cluster)
  colnames(result.sub)=c("cell","fraction","class")
  result=rbind(result,result.sub)
}
library(ggpubr)
plot.data=result
library(ggplot2)
library(ggprism)
library(ggsci)
plot.data$class=factor(plot.data$class,c("IM","EP"))
p=ggboxplot(plot.data, x="cell", y="fraction", color="class",
            xlab="",
            ylab="Fraction",
            title="Cibersort test",
            legend.title="Type",
            width=0.8,
            palette = c("#61b3de", "#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=class),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("Cibersort/compare_test.pdf",width=10,height=6)
p
dev.off()
##---immport test----
data=read.table('immport/immune_reaction_gsva_matrix.txt',sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
data=cbind(pathway=rownames(data),data)
data1=data %>%
  gather(sample, score, -pathway)
data=merge(data1,cluster,by="sample")
plot.data=data
plot.data$cluster=factor(plot.data$cluster,c("IM","EP"))
p=ggboxplot(plot.data, x="pathway", y="score", color="cluster",
            xlab="",
            ylab="score",
            title="ssGSEA test",
            legend.title="",
            width=0.8,
            palette = c("#61b3de", "#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=cluster),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("immport/test_compare.pdf",width=10,height=6)
p
dev.off()
###############
#----免疫细胞比较 整合在一起 只画显著的---------
rm(list=ls())
library(pheatmap)
library(colorspace)
library(plotrix)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(reshape2)
library(ggpubr)
##---免疫细胞比较 train------
setwd("D:/project/07.STAD_subtype/analysis0507_v7/08.immune_cell")
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
cluster=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", cluster$sample)
rownames(cluster)=cluster$sample
cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))

##--estimate
estimate=read.table("estimate/ESTIMATE_score.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names=F)
colnames(estimate)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(estimate))
sample=intersect(colnames(estimate),cluster$sample)
cluster=cluster[sample,]
cluster=cluster[order(cluster$cluster),]
estimate=estimate[-4,]
estimate=t(estimate)
estimate=estimate[cluster$sample,]
epic=read.table("epic/epic.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(epic)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(epic))
epic=dplyr::select(epic,c("cell_type",cluster$sample))
epic$cell_type=paste("epic:",epic$cell_type,sep="")
quantiseq=read.table("quantiseq/quantiseq.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(quantiseq)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(quantiseq))
quantiseq=dplyr::select(quantiseq,c("cell_type",cluster$sample))
quantiseq$cell_type=paste("quantiseq:",quantiseq$cell_type,sep="")
xcell=read.table("xcell/xcell.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(xcell)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(xcell))
xcell=dplyr::select(xcell,c("cell_type",cluster$sample))
xcell$cell_type=paste("xcell:",xcell$cell_type,sep="")
timer=read.table("timer/timer.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(timer)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(timer))
timer=dplyr::select(timer,c("cell_type",cluster$sample))
timer$cell_type=paste("timer:",timer$cell_type,sep="")
ssGSEA=read.table('ssGSEA/gsva_matrix.txt',sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(ssGSEA)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(ssGSEA))
ssGSEA=cbind(cell_type=rownames(ssGSEA),ssGSEA)
ssGSEA=dplyr::select(ssGSEA,c("cell_type",cluster$sample))
ssGSEA$cell_type=paste("ssGSEA:",ssGSEA$cell_type,sep="")
immune=rbind(epic,quantiseq,timer,xcell,ssGSEA)
data=melt(immune)
colnames(data)=c("cell_type","sample","value")
data=merge(data,cluster,by="sample")
K_W_test <- function(cell_type,data){
  loc=which(data$cell_type==cell_type)
  sub=data[loc,]
  aa=wilcox.test(value~cluster,data=sub)
  temp=data.frame(cell_type=cell_type,p=round(aa$p.value,4),z=aa$statistic)
  return(temp)
}
cell_type=unique(data$cell_type)
result=NULL
for( i in 1:length(cell_type)){
  temp=K_W_test(cell_type[i],data)
  result=rbind(result,temp)
}
sig=result[result$p<0.05,]
sig=sig[-grep("uncharacterized",sig$cell_type),]
sig=sig[-grep("score",sig$cell_type),]
immune=immune[immune$cell_type %in% sig$cell_type,]
rownames(immune)=immune$cell_type
plot.data=data[data$cell_type %in% sig$cell_type,]
p=ggboxplot(plot.data, x="cell_type", y="value", color="cluster",
            xlab="",
            ylab="Fraction",
            title="",
            legend.title="Type",
            width=0.8,
            palette = c("#61b3de","#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=cluster),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("immune_cell_compare.pdf",width=10,height=6)
p
dev.off()
cluster1=dplyr::select(immune,cluster$sample[cluster$cluster=="IM"])
cluster2=dplyr::select(immune,cluster$sample[cluster$cluster=="EP"])
cluster1_mean=apply(cluster1,1,sd)
cluster2_mean=apply(cluster2,1,sd)
aa=data.frame(IM=cluster1_mean,EP=cluster2_mean)
write.table(result,"wilcox_test_immune_cell.txt",sep="\t",quote=FALSE,row.names = F)
data1=immune[,-1]
pheno <- c(rep("IM", grep("IM",cluster$cluster)%>%length()),
           rep("EP", grep("EP",cluster$cluster)%>%length()))##样本的类型信息，为了后面画图的时候，加那个注释条
pheno=merge(cluster,estimate,by.x="sample",by.y="row.names")
pheno=pheno[order(pheno$cluster),]
#bk = unique(c(seq(0,1, length=50)))
normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}
data1=data1[,pheno$sample]
nor_data <- normalization(data1)
library(pheatmap)
bk <- c(seq(0,1,by=0.0001),seq(0,1,by=0.0001))
epic.length=length(grep("epic",rownames(nor_data)))
quantiseq.length=length(grep("quantiseq",rownames(nor_data)))
timer.length=length(grep("timer",rownames(nor_data)))
xcell.length=length(grep("xcell",rownames(nor_data)))
ssGSEA.length=length(grep("ssGSEA",rownames(nor_data)))
cluster1.length=length(which(cluster$cluster=="IM"))
annotation_row = data.frame(
  Method = c(rep("Epic",epic.length), rep("Quantiseq", quantiseq.length),
             rep("Timer",timer.length),rep("Xcell",xcell.length),
             rep("Zlatko.et.al",ssGSEA.length)))
row.names(annotation_row) <- rownames(data1)
annotation_col = pheno[,-1]
rownames(annotation_col)<-pheno$sample
StromalScore_col=color.scale(c(pheno$StromalScore,min(pheno$StromalScore),max(pheno$StromalScore)),na.color="#EDEDED55",extremes=c("#CCEBC5","#59A14F"))[1:length(pheno$StromalScore)]
ImmuneScore_col=color.scale(c(pheno$ImmuneScore,min(pheno$ImmuneScore),max(pheno$ImmuneScore)),na.color="#EDEDED55",extremes=c("#e4b1ab","#BC3C29FF"))[1:length(pheno$ImmuneScore)]
ESTIMATEScore_col=color.scale(c(pheno$ImmuneScore,min(pheno$ImmuneScore),max(pheno$ImmuneScore)),na.color="#EDEDED55",extremes=c("#FFFFB3","#F28E2B"))[1:length(pheno$ImmuneScore)]
colors=list(cluster=c(IM="#61b3de", EP="#ffee93"),
            StromalScore=StromalScore_col,
            ImmuneScore=ImmuneScore_col,
            ESTIMATEScore=ESTIMATEScore_col,
            Method=c("Epic"="#8DD3C7","Quantiseq"="#FFFFB3","Timer"="#BEBADA","Xcell"="#FB8072",
                     "Zlatko.et.al"="#FDB462"))
pdf(file="immune_cell_pheatmap.pdf",width=6,height =8)
pheatmap(data1,scale='row',
         show_colnames = F,#不展示行名
         cluster_rows = F,##对行不聚类
         cluster_cols = F,##对列不聚类
         annotation_col = annotation_col,##加注释
         annotation_row=annotation_row,
         annotation_colors = colors,
         gaps_col =cluster1.length,
         gaps_row = c(epic.length,epic.length+quantiseq.length,
                      epic.length+quantiseq.length+timer.length,
                      epic.length+quantiseq.length+timer.length+xcell.length
         ),
         color = colorRampPalette(c("navy","white","#BC3C29FF"))(100),
         #cellwidth=8,cellheight=8,##确定每个小格子的宽度和高度
         fontsize=8)##字体大小
dev.off()
library(circlize)
Method=annotation_row$Method
scaled_expr=t(scale(t(data1))) 
ha1=HeatmapAnnotation(df=pheno[,-1],
                      col=list(cluster=c(IM="#61b3de", EP="#ffee93"),
                               StromalScore= colorRamp2(c(min(pheno$StromalScore),max(pheno$StromalScore)), c("#CCEBC5","#59A14F")),
                               ImmuneScore= colorRamp2(c(min(pheno$ImmuneScore),max(pheno$ImmuneScore)), c("#e4b1ab","#BC3C29FF")),
                               ESTIMATEScore=colorRamp2(c(min(pheno$ESTIMATEScore),max(pheno$ESTIMATEScore)), c("#FFFFB3","#F28E2B"))))
ha2=rowAnnotation(df=data.frame(
  Method = c(rep("Epic",epic.length), rep("Quantiseq", quantiseq.length),
             rep("Timer",timer.length),rep("Xcell",xcell.length),
             rep("Zlatko.et.al",ssGSEA.length))),
  col=list(Method=c("Epic"="#8DD3C7","Quantiseq"="#FFFFB3","Timer"="#BEBADA","Xcell"="#FB8072",
                    "Zlatko.et.al"="#FDB462")))


annotation_row = data.frame(
  Method = c(rep("Epic",epic.length), rep("Quantiseq", quantiseq.length),
             rep("Timer",timer.length),rep("Xcell",xcell.length),
             rep("Zlatko.et.al",ssGSEA.length)))
row_labels=str_split_fixed(rownames(scaled_expr),":",2)[,2]
pdf(file="immune_cell_complexheatmap_nocluster_1.pdf",width=8,height =8)
Heatmap(scaled_expr,name="z-score",top_annotation=ha1,left_annotation =ha2 ,
        cluster_rows=TRUE,cluster_columns=FALSE,show_row_names=TRUE,
        show_column_names=FALSE,
        column_split = cluster$cluster,
        row_split = Method,
        row_labels = row_labels,
        #column_names_gp=gpar(fontsize=5),
        row_names_gp = gpar(fontsize=8))
dev.off()
##---免疫细胞比较 test------
rm(list=ls())
library(pheatmap)
library(colorspace)
library(plotrix)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(reshape2)
library(ggpubr)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/08.immune_cell")
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
cluster=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", cluster$sample)
rownames(cluster)=cluster$sample
cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
##--estimate
estimate=read.table("estimate/ESTIMATE_score.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names=F)
colnames(estimate)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(estimate))
sample=intersect(colnames(estimate),cluster$sample)
cluster=cluster[sample,]
cluster=cluster[order(cluster$cluster),]
estimate=estimate[-4,]
estimate=t(estimate)
estimate=estimate[cluster$sample,]
epic=read.table("epic/epic.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(epic)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(epic))
epic=dplyr::select(epic,c("cell_type",cluster$sample))
epic$cell_type=paste("epic:",epic$cell_type,sep="")
quantiseq=read.table("quantiseq/quantiseq.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(quantiseq)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(quantiseq))
quantiseq=dplyr::select(quantiseq,c("cell_type",cluster$sample))
quantiseq$cell_type=paste("quantiseq:",quantiseq$cell_type,sep="")
xcell=read.table("xcell/xcell.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(xcell)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(xcell))
xcell=dplyr::select(xcell,c("cell_type",cluster$sample))
xcell$cell_type=paste("xcell:",xcell$cell_type,sep="")
timer=read.table("timer/timer.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(timer)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(timer))
timer=dplyr::select(timer,c("cell_type",cluster$sample))
timer$cell_type=paste("timer:",timer$cell_type,sep="")
ssGSEA=read.table('ssGSEA/gsva_matrix.txt',sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(ssGSEA)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(ssGSEA))
ssGSEA=cbind(cell_type=rownames(ssGSEA),ssGSEA)
ssGSEA=dplyr::select(ssGSEA,c("cell_type",cluster$sample))
ssGSEA$cell_type=paste("ssGSEA:",ssGSEA$cell_type,sep="")
immune=rbind(epic,quantiseq,timer,xcell,ssGSEA)
data=melt(immune)
colnames(data)=c("cell_type","sample","value")
data=merge(data,cluster,by="sample")
K_W_test <- function(cell_type,data){
  loc=which(data$cell_type==cell_type)
  sub=data[loc,]
  aa=wilcox.test(value~cluster,data=sub)
  tEPp=data.frame(cell_type=cell_type,p=round(aa$p.value,4),z=aa$statistic)
  return(tEPp)
}
cell_type=unique(data$cell_type)
result=NULL
for( i in 1:length(cell_type)){
  tEPp=K_W_test(cell_type[i],data)
  result=rbind(result,tEPp)
}
sig=result[result$p<0.05,]
sig=sig[-grep("uncharacterized",sig$cell_type),]
sig=sig[-grep("score",sig$cell_type),]
immune=immune[immune$cell_type %in% sig$cell_type,]
rownames(immune)=immune$cell_type
plot.data=data[data$cell_type %in% sig$cell_type,]
p=ggboxplot(plot.data, x="cell_type", y="value", color="cluster",
            xlab="",
            ylab="Fraction",
            title="",
            legend.title="Type",
            width=0.8,
            palette = c("#61b3de","#ffee93"))+
  rotate_x_text(90)+ ##这个是改倾斜度的
  stat_compare_means(aes(group=cluster),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf("immune_cell_compare_test.pdf",width=10,height=6)
p
dev.off()
cluster1=dplyr::select(immune,cluster$sample[cluster$cluster=="IM"])
cluster2=dplyr::select(immune,cluster$sample[cluster$cluster=="EP"])
cluster1_mean=apply(cluster1,1,sd)
cluster2_mean=apply(cluster2,1,sd)
aa=data.frame(IM=cluster1_mean,EP=cluster2_mean)
write.table(result,"wilcox_test_immune_cell.txt",sep="\t",quote=FALSE,row.names = F)
data1=immune[,-1]
pheno <- c(rep("IM", grep("IM",cluster$cluster)%>%length()),
           rep("EP", grep("EP",cluster$cluster)%>%length()))##样本的类型信息，为了后面画图的时候，加那个注释条
pheno=merge(cluster,estimate,by.x="sample",by.y="row.names")
pheno=pheno[order(pheno$cluster),]
#bk = unique(c(seq(0,1, length=50)))
normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}
data1=data1[,pheno$sample]
nor_data <- normalization(data1)
library(pheatmap)
bk <- c(seq(0,1,by=0.0001),seq(0,1,by=0.0001))
epic.length=length(grep("epic",rownames(nor_data)))
quantiseq.length=length(grep("quantiseq",rownames(nor_data)))
timer.length=length(grep("timer",rownames(nor_data)))
xcell.length=length(grep("xcell",rownames(nor_data)))
ssGSEA.length=length(grep("ssGSEA",rownames(nor_data)))
cluster1.length=length(which(cluster$cluster=="IM"))
annotation_row = data.frame(
  Method = c(rep("Epic",epic.length), rep("Quantiseq", quantiseq.length),
             rep("Timer",timer.length),rep("Xcell",xcell.length),
             rep("Zlatko.et.al",ssGSEA.length)))
row.names(annotation_row) <- rownames(data1)
annotation_col = pheno[,-1]
rownames(annotation_col)<-pheno$sample
StromalScore_col=color.scale(c(pheno$StromalScore,min(pheno$StromalScore),max(pheno$StromalScore)),na.color="#EDEDED55",extremes=c("#CCEBC5","#59A14F"))[1:length(pheno$StromalScore)]
ImmuneScore_col=color.scale(c(pheno$ImmuneScore,min(pheno$ImmuneScore),max(pheno$ImmuneScore)),na.color="#EDEDED55",extremes=c("#e4b1ab","#BC3C29FF"))[1:length(pheno$ImmuneScore)]
ESTIMATEScore_col=color.scale(c(pheno$ImmuneScore,min(pheno$ImmuneScore),max(pheno$ImmuneScore)),na.color="#EDEDED55",extremes=c("#FFFFB3","#F28E2B"))[1:length(pheno$ImmuneScore)]
colors=list(cluster=c(IM="#61b3de", EP="#ffee93"),
            StromalScore=StromalScore_col,
            ImmuneScore=ImmuneScore_col,
            ESTIMATEScore=ESTIMATEScore_col,
            Method=c("Epic"="#8DD3C7","Quantiseq"="#FFFFB3","Timer"="#BEBADA","Xcell"="#FB8072",
                     "Zlatko.et.al"="#FDB462"))
pdf(file="immune_cell_pheatmap_test.pdf",width=6,height = 6)
pheatmap(data1,scale='row',
         show_colnames = F,#不展示行名
         cluster_rows = F,##对行不聚类
         cluster_cols = F,##对列不聚类
         annotation_col = annotation_col,##加注释
         annotation_row=annotation_row,
         annotation_colors = colors,
         gaps_col =cluster1.length,
         gaps_row = c(epic.length,epic.length+quantiseq.length,
                      epic.length+quantiseq.length+timer.length,
                      epic.length+quantiseq.length+timer.length+xcell.length
         ),
         color = colorRampPalette(c("navy","white","#BC3C29FF"))(100),
         #cellwidth=8,cellheight=8,##确定每个小格子的宽度和高度
         fontsize=8)##字体大小
dev.off()
library(circlize)
Method=annotation_row$Method
scaled_expr=t(scale(t(data1))) 
ha1=HeatmapAnnotation(df=pheno[,-1],
                      col=list(cluster=c(IM="#61b3de", EP="#ffee93"),
                               StromalScore= colorRamp2(c(min(pheno$StromalScore),max(pheno$StromalScore)), c("#CCEBC5","#59A14F")),
                               ImmuneScore= colorRamp2(c(min(pheno$ImmuneScore),max(pheno$ImmuneScore)), c("#e4b1ab","#BC3C29FF")),
                               ESTIMATEScore=colorRamp2(c(min(pheno$ESTIMATEScore),max(pheno$ESTIMATEScore)), c("#FFFFB3","#F28E2B"))))
ha2=rowAnnotation(df=data.frame(
  Method = c(rep("Epic",epic.length), rep("Quantiseq", quantiseq.length),
             rep("Timer",timer.length),rep("Xcell",xcell.length),
             rep("Zlatko.et.al",ssGSEA.length))),
  col=list(Method=c("Epic"="#8DD3C7","Quantiseq"="#FFFFB3","Timer"="#BEBADA","Xcell"="#FB8072",
                    "Zlatko.et.al"="#FDB462")))


annotation_row = data.frame(
  Method = c(rep("Epic",epic.length), rep("Quantiseq", quantiseq.length),
             rep("Timer",timer.length),rep("Xcell",xcell.length),
             rep("Zlatko.et.al",ssGSEA.length)))
row_labels=str_split_fixed(rownames(scaled_expr),":",2)[,2]
pdf(file="immune_cell_complexheatmap_nocluster_test_1.pdf",width=8,height =8)
Heatmap(scaled_expr,name="z-score",top_annotation=ha1,left_annotation =ha2 ,
        cluster_rows=TRUE,cluster_columns=FALSE,show_row_names=TRUE,
        show_column_names=FALSE,
        column_split = cluster$cluster,
        row_split = Method,
        row_labels = row_labels,#自定义行名标签
        #column_names_gp=gpar(fontsize=5),
        row_names_gp = gpar(fontsize=8))
dev.off()
#############
##---临床特征和亚型的关系------
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(stringr)
library(cowplot)
###--train-------------
clinical_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/clinical"
clinical=read.table(file.path(clinical_path,"tumor_clinical.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
cluster=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", cluster$sample)
rownames(cluster)=cluster$sample
cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
data=merge(clinical,cluster,by="sample")
data1=data
clustercolor <- c("#61b3de", "#ffee93") 
names(clustercolor) <- c("IM","EP") #类型颜色
Agecolor <- color.scale(c(clinical$age,35,85),na.color="#EDEDED55",extremes=c("#e4b1ab","#cc444b"))[1:length(clinical$age)]
names(Agecolor)<- clinical$age
gendercolor <- c("#ffee93", "#a0ced9") 
names(gendercolor) <- c("female","male") #类型颜色
Tcolor <- c("#ebebeb", "#caf0f8",  "#00b4d8","#277da1")
names(Tcolor) <- c("T1","T2","T3","T4") #类型颜色
Mcolor <- c("#8dd3c7","#ffffb3")
names(Mcolor) <- c("M0","M1") #类型颜色
Ncolor <- c("#ffffb3", "#42b5407f","#0099B47F","#00468B7F" )
names(Ncolor) <- c("N0","N1","N2","N3") #类型颜色
stagecolor <- c("#bfd7ea", "#4E80AB","#a7c957", "#6a994e")
names(stagecolor) <- c("stage i","stage ii","stage iii","stage iv") #类型颜色
gradecolor <- c("#bbbdc0", "#f6ddb4","#f6a34a")
names(gradecolor) <- c("G1","G2","G3") #类型颜色
ann_colors <- list(Cluster=clustercolor, Age= Agecolor, Gender=gendercolor, Stage=stagecolor,Grade=gradecolor,
                   T=Tcolor,N=Ncolor,M=Mcolor) #颜色设置
##--分类的临床特征和亚型的关系--#
fisher_plot <- function(data1,character){
  ##---计算fisher检验得到的p值--#
  data1=data1[!is.na(data1[,character]),c('cluster',character)]
  table_clinical=table(data1[,character],data1[,'cluster'])
  fisher_res=fisher.test(table_clinical)
  pv=round(fisher_res$p.value,3)
  colnames(data1)[2]="clinical"
  character1=str_to_title(character)
  data1$cluster=factor(data1$cluster,levels=c("IM","EP"))
  if(pv<0.001){pv="P<0.001"}else{pv=paste("P=",round(pv,3),sep="")}
  color=ann_colors[[character1]]
  data1$clinical=factor(data1$clinical,levels=names(color))
  IM_num=length(which(data1$cluster=="IM"))
  EM_num=length(which(data1$cluster=="EP"))
  p=ggplot(data1,aes(x=cluster))+ geom_bar(aes(fill=clinical),width=0.8, position="fill")+
    labs(y="Percentage of Patients",x="")+scale_fill_manual(values = color)+
    annotate(x=1.5,y=1.05,geom = "text",label=pv,size=5)+
    theme_bw()+
    #theme(axis.text.x = element_text(angle = 0,vjust = 0.85,hjust = 0.75))+
    guides(fill = guide_legend(title = character1))+
    scale_x_discrete(labels = c(paste0("IM","\n","n=(",IM_num,")"),paste0("EP","\n","n=(",EM_num,")")))+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.y = element_text(size = 12,
                                     face="bold"))
  
  pdf(paste0(character,"_barplot.pdf"),height=5,width=4)
  print(p)
  dev.off()
  return(p)
}
##--连续的临床特征和亚型的关系--#
boxplot_plot <- function(data1,character){
  data1=data1[!is.na(data1[,character]),c('cluster',character)]
  colnames(data1)[2]="clinical"
  character1=str_to_title(character)
  data1$cluster=factor(data1$cluster,levels=c("IM","EP"))
  IM_num=length(which(data1$cluster=="IM"))
  EM_num=length(which(data1$cluster=="EP"))
  p=ggplot(data1,aes(x=cluster,y=clinical,fill=cluster))+geom_boxplot(width=0.8)+
    labs(y=character1,x="")+scale_fill_manual(values = c("#61b3de", "#ffee93"))+
    stat_compare_means(aes(group = cluster), 
                       method = "wilcox.test",label.x=1.25)+
    theme_bw()+
    #theme(axis.text.x = element_text(angle = 0,vjust = 0.85,hjust = 0.75))+
    theme(legend.position = 'none')+
    scale_x_discrete(labels = c(paste0("IM","\n","n=(",IM_num,")"),paste0("EP","\n","n=(",EM_num,")")))+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.y = element_text(size = 12,
                                     face="bold"))
  
  pdf(paste0(character,"_boxplot.pdf"),height=5,width=3.5)
  print(p)
  dev.off()
  return(p)
}
setwd("D:/project/07.STAD_subtype/analysis0507_v7//10.clinical")
character=colnames(data)[-c(1,9)]
p_list=list()
for(i in 2:length(character)){
  # if(character[i]=="age"){
  #   p_list[[i]]=boxplot_plot(data,character[i])
  # }else{
  #   p_list[[i]]=fisher_plot(data,character[i])
  # }
  p_list[[i-1]]=fisher_plot(data,character[i])
}
boxplot_plot(data,"age")
pdf("clinical_cluster.pdf",width=6,height=12)
plot_grid(plotlist=p_list,ncol=2)
dev.off()
#########
###--test-------------
clinical_path="D:/project/07.STAD_subtype/analysis0507_v7//00.data/clinical"
clinical=read.table(file.path(clinical_path,"tumor_clinical.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
cluster=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", cluster$sample)
rownames(cluster)=cluster$sample
cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
data=merge(clinical,cluster,by="sample")
data1=data
setwd("D:/project/07.STAD_subtype/analysis0507_v7//10.clinical")
dir.create("test")
setwd("test")
character=colnames(data)[-c(1,9)]
p_list_test=list()
for(i in 2:length(character)){
  # if(character[i]=="age"){
  #   p_list[[i]]=boxplot_plot(data,character[i])
  # }else{
  #   p_list[[i]]=fisher_plot(data,character[i])
  # }
  p_list_test[[i-1]]=fisher_plot(data,character[i])
}
boxplot_plot(data,"age")
pdf("clinical_cluster.pdf",width=6,height=12)
plot_grid(plotlist=p_list_test,ncol=2)
dev.off()

##---年龄和亚型的关系----##
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(stringr)
library(cowplot)
###--train-------------
clinical_path="D:/project/07.zhouY_project/01.analysis0422/00.data/clinical"
clinical=read.table(file.path(clinical_path,"tumor_clinical.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster_path="D:/project/07.zhouY_project/01.analysis0422/04.important_TF_cox_CpG_cor_again/GSVA_methylation/"
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
train=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
train$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", train$sample)
rownames(train)=train$sample
train$cluster=factor(train$cluster,levels=c("IM","EP"))
test$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", test$sample)
rownames(test)=test$sample
test$cluster=factor(test$cluster,levels=c("IM","EP"))
train$dataset="train"
test$dataset="test"
cluster=rbind(train,test)
data=merge(clinical,cluster,by="sample")
data1=data
data1=data1[!is.na(data1[,"age"]),]
data_long=data1[,c("sample",'cluster',"age","dataset")]
colnames(data_long)=c("sample","group2","values","cancer")
data_long$group=paste0(data_long$group2,data_long$cancer)
data_long$group=factor(data_long$group,levels=c("IMtrain","EPtrain","IMtest","EPtest"))
data1$cluster=factor(data1$cluster,levels=c("IM","EP"))
data1$dataset=factor(data1$dataset,levels=c("train","test"))
number=as.data.frame(table(data1[,c("cluster","dataset")]))
p <- ggplot(data_long)+
  geom_boxplot(aes(group, values, fill = group2), 
               # 间距调整：
               width = 0.5,
               outlier.shape = NA)+
  labs(y="Age",x="")+
  # 顶部灰色方块：
  geom_rect(aes(xmin=0, xmax=5, ymin=95, ymax = 100), fill="#eaeae0")+
  # 灰色竖线：
  geom_vline(xintercept = c(2.5), color = "#bcbdbf", alpha = 0.8)+
  scale_x_discrete(labels = c(paste0("IM","\n","n=(",number$Freq[1],")"),
                              paste0("EP","\n","n=(",number$Freq[2],")"),
                              paste0("IM","\n","n=(",number$Freq[3],")"),
                              paste0("EP","\n","n=(",number$Freq[4],")")
  ))+
  scale_y_continuous(expand = c(0,0))+
  # 颜色：
  scale_fill_manual(values = c("#61b3de", "#ffee93"))+
  # 主题：
  theme_bw()+
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))+
  # 添加p值：
  geom_signif(aes(group, values),
              comparisons = list(c("IMtrain", "EPtrain"),
                                 c("IMtest", "EPtest")),
              vjust = 2, 
              tip_length = rep(0,2), 
              y_position = c(90, 90),
              # 修改注释内容
              annotations = c("Wilcoxon,p = 0.023", "Wilcoxon,p = 0.015"))+
  annotate("text", x = c(1.5,3.5), y = 98, label = c("Train Datesets","Test Datesets"))
setwd("D:/project/07.STAD_subtype/analysis0507_v7/10.clinical")
pdf("age_boxplot_v2.pdf",width=4,height=4.5)
p
dev.off()

#####################

##---年龄和亚型的关系----
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(stringr)
library(cowplot)
clinical_path="D:/project/07.zhouY_project/01.analysis0422/00.data/clinical"
clinical=read.table(file.path(clinical_path,"tumor_clinical.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
train=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
train$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", train$sample)
rownames(train)=train$sample
train$cluster=factor(train$cluster,levels=c("IM","EP"))
test$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", test$sample)
rownames(test)=test$sample
test$cluster=factor(test$cluster,levels=c("IM","EP"))
train$dataset="Train Datasets"
test$dataset="Test Datasets"
cluster=rbind(train,test)
data=merge(clinical,cluster,by="sample")
data1=data
data1=data1[!is.na(data1[,"age"]),]
data1$cluster=factor(data1$cluster,levels=c("IM","EP"))
data1$dataset=factor(data1$dataset,levels=c("Train Datasets","Test Datasets"))
data_long=data1[,c("sample",'cluster',"age","dataset")]
colnames(data_long)=c("sample","group2","values","cancer")
data_long$group=paste0(data_long$group2,data_long$cancer)
data_long$group=factor(data_long$group,levels=c("IMtrain","EPtrain","IMtest","EPtest"))
number=as.data.frame(table(data1[,c("cluster","dataset")]))
p=ggboxplot(data_long, x="group2", y="values", color="group2",
            xlab="",
            ylab="Age",
            legend.title="Type",
            width=0.5,
            facet.by="cancer", #分面
            palette = c("#61b3de", "#ffee93"),
            add = "jitter")+
  rotate_x_text(0)+ ##这个是改倾斜度的
  
  stat_compare_means(aes(group = group2), 
                     method = "wilcox.test")+
  scale_x_discrete(labels = c(paste0("IM","\n","n=(",number$Freq[1],")"),
                              paste0("EP","\n","n=(",number$Freq[2],")"),
                              paste0("IM","\n","n=(",number$Freq[3],")"),
                              paste0("EP","\n","n=(",number$Freq[4],")")
  ))
setwd("D:/project/07.STAD_subtype/analysis0507_v7/10.clinical")
pdf("age_boxplot_v3.pdf",width=4,height=4.5)
p
dev.off()


##---MSI 亚组比较-------
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(stringr)
library(cowplot)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/11.MSI")
MSI=read.table("stad_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",header=T,stringsAsFactors = F,quote="")
MSI=MSI[,c("Patient.ID","MSI.MANTIS.Score","MSIsensor.Score")]
colnames(MSI)=c("sample","MSI_socre1","MSI_score2")
MSI=MSI[!is.na(MSI$MSI_socre1),]
#MSI Score reported by MANTIS. The suggested thresholds are MSI: >0.6(MSI-H), Indeterminate: 0.4-0.6(MSS) and MSS: <0.4(MSI-L).
#MSI Score reported by MSIsensor. The suggested thresholds are MSI: >10, Indeterminate: 4-10 and MSS: <10.
MSI$MANTIS="MSS"
MSI$MANTIS[MSI$MSI_socre1>0.6]="MSI-H"
MSI$MANTIS[MSI$MSI_socre1<0.4]="MSI-L"
MSI$MSIsensor="MSS"
MSI$MSIsensor[MSI$MSI_score2>10]="MSI-H"
MSI$MSIsensor[MSI$MSI_score2<4]="MSI-L"
write.table(MSI,"MSI.xls",sep="\t",quote=FALSE,row.names = F)
clinical_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/clinical"
clinical=read.table(file.path(clinical_path,"tumor_clinical.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
train=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
train$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", train$sample)
rownames(train)=train$sample
train$cluster=factor(train$cluster,levels=c("IM","EP"))
test$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", test$sample)
rownames(test)=test$sample
test$cluster=factor(test$cluster,levels=c("IM","EP"))
clinical_barplot <- function(cluster,MSI,character,title){
  data=merge(cluster,MSI,by="sample")
  data1=data[,c("cluster",character)]
  colnames(data1)[2]="clinical"
  table_clinical=table(data1[,"clinical"],data1[,'cluster'])
  fisher_res=fisher.test(table_clinical)
  pv=round(fisher_res$p.value,3)
  if(pv<0.001){pv="P<0.001"}else{pv=paste("P=",round(pv,3),sep="")}
  color=c("#925fa7", "#c5b3d1", "#bbbdc0")
  data1$clinical=factor(data1$clinical,levels=c("MSI-H","MSI-L","MSS"))
  cluster1_num=length(which(data1$cluster=="IM"))
  cluster2_num=length(which(data1$cluster=="EP"))
  p=ggplot(data1,aes(x=cluster))+ geom_bar(aes(fill=clinical),width=0.8, position="fill")+
    labs(y="Percentage of Patients",x="",title=title)+scale_fill_manual(values = color)+
    annotate(x=1.5,y=1.05,geom = "text",label=pv,size=5)+
    theme_bw()+
    #theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))+
    guides(fill = guide_legend(title = character))+
    scale_x_discrete(labels = c(paste0("IM","\n","n=(",cluster1_num,")"),paste0("EP","\n","n=(",cluster2_num,")")))
  
  return(p)
}

p1_train=clinical_barplot(train,MSI,"MANTIS","Train Datasets")
p2_train=clinical_barplot(train,MSI,"MSIsensor","Train Datasets")
p_list_train=list(p1_train,p2_train)
p1_test=clinical_barplot(test,MSI,"MANTIS","Test Datasets")
p2_test=clinical_barplot(test,MSI,"MSIsensor","Test Datasets")
p_list_test=list(p1_test,p2_test)
p_list1=list(p1_train,p1_test)
p_list2=list(p2_train,p2_test)
pdf("train_MSI.pdf",height=4.5,width=6)
plot_grid(plotlist=p_list_train,ncol=2)
dev.off()
pdf("test_MSI.pdf",height=4.5,width=6)
plot_grid(plotlist=p_list_test,ncol=2)
dev.off()
pdf("MANTIS_MSI.pdf",height=4.5,width=6)
plot_grid(plotlist=p_list1,ncol=2)
dev.off()
pdf("MSIsensor_MSI.pdf",height=4.5,width=6)
plot_grid(plotlist=p_list2,ncol=2)
dev.off()
##---TCGA 已有分型和现有分型比-------
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(stringr)
library(cowplot)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/11.MSI")
clinical=read.table("stad_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",header=T,stringsAsFactors = F,quote="")
clinical=clinical[,c("Patient.ID","Subtype")]
colnames(clinical)=c("sample","Subtype")
clinical=clinical[clinical$Subtype!="STAD_POLE",]
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
train=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
train$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", train$sample)
rownames(train)=train$sample
train$cluster=factor(train$cluster,levels=c("IM","EP"))
test$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", test$sample)
rownames(test)=test$sample
test$cluster=factor(test$cluster,levels=c("IM","EP"))
cluster=rbind(train,test)
clinical_barplot <- function(cluster,clinical,title){
  data=merge(cluster,clinical,by="sample")
  data1=data[,c("cluster","Subtype")]
  data1=data1[!is.na(data1$Subtype),]
  table_clinical=table(data1[,'cluster'],data1[,"Subtype"])
  fisher_res=fisher.test(table_clinical)
  pv=round(fisher_res$p.value,3)
  data1$cluster=factor(data1$cluster,levels=c("IM","EP"))
  if(pv<0.001){pv="P<0.001"}else{pv=paste("P=",round(pv,3),sep="")}
  p=ggplot(data1,aes(x=Subtype))+ geom_bar(aes(fill=cluster),width=0.8, position="fill")+
    labs(y="Percentage of Patients",x="",title=title)+
    scale_fill_manual(values = c("#61b3de","#ffee93"))+
    annotate(x=3,y=1.05,geom = "text",label=pv,size=5)+
    theme_bw()+theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))
  
  return(p)
}

p1_train=clinical_barplot(train,clinical,"Train Datasets")
p1_test=clinical_barplot(test,clinical,"Test Datasets")
p3=clinical_barplot(cluster,clinical,"All Datasets")
p_list1=list(p1_train,p1_test,p3)
pdf("subtype_type.pdf",height=4.5,width=12)
plot_grid(plotlist=p_list1,ncol=3)
dev.off()
############################
##---clinical pheatmap add MSI-----
rm(list=ls())
library(pheatmap)
library(colorspace)
library(plotrix)
library(RColorBrewer)
library(stringr)
data_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation"
train_score=read.table(file.path(data_path,"train_TF_cg_score.xls"),header=T,sep="\t",stringsAsFactors = F,row.names = 1,check.names = F)##
test_score=read.table(file.path(data_path,"test_TF_cg_score.xls"),header=T,sep="\t",stringsAsFactors = F,row.names = 1,check.names = F)##
colnames(train_score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(train_score))
colnames(test_score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(test_score))
TF_path="D:/project/07.STAD_subtype/analysis0507_v7/04.OS/TF"
TF=read.table(file.path(TF_path,"single_cox_sig.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
train=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
train$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", train$sample)
rownames(train)=train$sample
train$cluster=factor(train$cluster,levels=c("IM","EP"))
test$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", test$sample)
rownames(test)=test$sample
test$cluster=factor(test$cluster,levels=c("IM","EP"))
train_cluster=train
test_cluster=test
clinical_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/clinical"
clinical=read.table(file.path(clinical_path,"tumor_clinical.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
MSI_path="D:/project/07.STAD_subtype/analysis0507_v7/11.MSI"
MSI=read.table(file.path(MSI_path,"MSI.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
MSI=MSI[,c("sample","MSIsensor")]
colnames(MSI)[2]="MSI"
clinical=merge(clinical,MSI,by="sample")
network_path="D:/project/07.STAD_subtype/analysis0507_v7/05.cox_TF_cox_CpG_cor/network_regulon"
network=read.table(file.path(network_path,"TF_sig_cox_regulon_cor.txt"),sep='\t',header=T,stringsAsFactors = F,quote="")
network=network[order(network$number,network$abs,decreasing = T),]
network1=network[,c("TF","regulon")]
pheatmap_clincial <- function(score,TF,clinical,cluster,network1){
  data=score[TF$Characteristics,]
  clinical=merge(clinical,cluster,by="sample")
  clinical$cluster=factor(clinical$cluster,levels=c("IM","EP"))
  clinical=clinical[order(clinical$cluster),]
  data=data[,clinical$sample]
  data=merge(network1,data,by.x="TF",by.y="row.names")
  rownames(data)=data$regulon
  data=data[,-c(1,2)]
  data=data[network1$regulon,]
  ##---color---##
  clustercolor <- c( "#61b3de","#ffee93") 
  names(clustercolor) <- c("IM","EP") #类型颜色
  Agecolor <- color.scale(c(clinical$age,35,85),na.color="#EDEDED55",extremes=c("#e4b1ab","#cc444b"))[1:length(clinical$age)]
  names(Agecolor)<- clinical$age
  gendercolor <- c("#ffee93", "#a0ced9") 
  names(gendercolor) <- c("female","male") #类型颜色
  Tcolor <- c("#ebebeb", "#caf0f8",  "#00b4d8","#277da1")
  names(Tcolor) <- c("T1","T2","T3","T4") #类型颜色
  Mcolor <- c("#8dd3c7","#ffffb3", "white")
  names(Mcolor) <- c("M0","M1",NA) #类型颜色
  Ncolor <- c("#ffffb3", "#42b5407f","#0099B47F","#00468B7F", "white")
  names(Ncolor) <- c("N0","N1","N2","N3",NA) #类型颜色
  stagecolor <- c("#bfd7ea", "#4E80AB","#a7c957", "#6a994e","white")
  names(stagecolor) <- c("stage i","stage ii","stage iii","stage iv",NA) #类型颜色
  gradecolor <- c("#bbbdc0", "#f6ddb4","#f6a34a", "white")
  names(gradecolor) <- c("G1","G2","G3",NA) #类型颜色
  MSIcolor <- c("#925fa7", "#c5b3d1", "#bbbdc0")
  names(MSIcolor) <- c("MSI-H","MSI-L","MSS") #类型颜色
  ann_colors <- list(Cluster=clustercolor, Age= Agecolor, MSI=MSIcolor,Gender=gendercolor, Stage=stagecolor,Grade=gradecolor,
                     T=Tcolor,N=Ncolor,M=Mcolor) #颜色设置
  # ann_colors <- list("Cluster"=clustercolor, 
  #                    Age= Agecolor, 
  #                    "Gender"=gendercolor, "Stage"=stagecolor,"Grade"=gradecolor,
  #                    "T"=Tcolor,"N"=Ncolor,"M"=Mcolor) #颜色设置
  colnames(clinical)=str_to_title(colnames(clinical))
  colnames(clinical)[which(colnames(clinical)=="Msi")]="MSI"
  annotation_col = clinical[,c("Cluster","Age","MSI","Gender","T","N","M","Stage","Grade")]
  # annotation_col = clinical[,c("Cluster",
  #                              "Age",
  #                              "Gender","T","N","M","Stage","Grade")]
  
  row.names(annotation_col) <- clinical$Sample
  #在画图的时候，heatmap函数中多添加annotation_colors即可：
  p=pheatmap(as.matrix(data),cluster_rows = T,cluster_cols = T,
             color=colorRampPalette(c("navy","white","firebrick3"))(100),
             show_colnames = F,border_color = NA,scale = "row",show_rownames =T,
             annotation_col = annotation_col,# annotation_row = annotation_row,
             annotation_colors = ann_colors)
  p1=pheatmap(data,cluster_rows =T,cluster_cols = F,
              color=colorRampPalette(c("navy","white","firebrick3"))(100),
              show_colnames = F,border_color = NA,scale = "row",show_rownames =T,
              annotation_col = annotation_col,# annotation_row = annotation_row,
              annotation_colors = ann_colors)
  p_list <- list(p,p1)
  return(p_list)
}
p_list=pheatmap_clincial(train_score,TF,clinical,train_cluster,network1)
p_list_test=pheatmap_clincial(test_score,TF,clinical,test_cluster,network1)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/10.clinical")
pdf("clinical_data_MSI_pheatmap.pdf",width=8,height=6.5)
p_list[[1]]
dev.off()
pdf("clinical_data_MSI_pheatmap_no_cluster.pdf",width=8,height=6.5)
p_list[[2]]
dev.off()
pdf("test/clinical_data_MSI_pheatmap.pdf",width=8,height=6.5)
p_list_test[[1]]
dev.off()
pdf("test/clinical_data_MSI_pheatmap_no_cluster.pdf",width=8,height=6.5)
p_list_test[[2]]
dev.off()
##--------no MSI---------
rm(list=ls())
##---clinical pheatmap---#
library(pheatmap)
library(colorspace)
library(plotrix)
library(RColorBrewer)
library(stringr)
data_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation"
train_score=read.table(file.path(data_path,"train_TF_cg_score.xls"),header=T,sep="\t",stringsAsFactors = F,row.names = 1,check.names = F)##
test_score=read.table(file.path(data_path,"test_TF_cg_score.xls"),header=T,sep="\t",stringsAsFactors = F,row.names = 1,check.names = F)##
colnames(train_score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(train_score))
colnames(test_score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(test_score))
TF_path="D:/project/07.STAD_subtype/analysis0507_v7/04.OS/TF"
TF=read.table(file.path(TF_path,"single_cox_sig.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
train=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
train$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", train$sample)
rownames(train)=train$sample
train$cluster=factor(train$cluster,levels=c("IM","EP"))
test$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", test$sample)
rownames(test)=test$sample
test$cluster=factor(test$cluster,levels=c("IM","EP"))
train_cluster=train
test_cluster=test
clinical_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/clinical"
clinical=read.table(file.path(clinical_path,"tumor_clinical.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
# MSI_path="D:/project/07.STAD_subtype/analysis0507_v7/11.MSI"
# MSI=read.table(file.path(MSI_path,"MSI.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
# MSI=MSI[,c("sample","MSIsensor")]
# colnames(MSI)[2]="MSI"
# clinical=merge(clinical,MSI,by="sample")
network_path="D:/project/07.STAD_subtype/analysis0507_v7/05.cox_TF_cox_CpG_cor/network_regulon"
network=read.table(file.path(network_path,"TF_sig_cox_regulon_cor.txt"),sep='\t',header=T,stringsAsFactors = F,quote="")
network=network[order(network$number,network$abs,decreasing = T),]
network1=network[,c("TF","regulon")]
pheatmap_clincial <- function(score,TF,clinical,cluster,network1){
  data=score[TF$Characteristics,]
  clinical=merge(clinical,cluster,by="sample")
  clinical$cluster=factor(clinical$cluster,levels=c("IM","EP"))
  clinical=clinical[order(clinical$cluster),]
  data=data[,clinical$sample]
  data=merge(network1,data,by.x="TF",by.y="row.names")
  rownames(data)=data$regulon
  data=data[,-c(1,2)]
  data=data[network1$regulon,]
  ##---color---##
  clustercolor <- c( "#61b3de","#ffee93") 
  names(clustercolor) <- c("IM","EP") #类型颜色
  Agecolor <- color.scale(c(clinical$age,35,85),na.color="#EDEDED55",extremes=c("#e4b1ab","#cc444b"))[1:length(clinical$age)]
  names(Agecolor)<- clinical$age
  gendercolor <- c("#ffee93", "#a0ced9") 
  names(gendercolor) <- c("female","male") #类型颜色
  Tcolor <- c("#ebebeb", "#caf0f8",  "#00b4d8","#277da1")
  names(Tcolor) <- c("T1","T2","T3","T4") #类型颜色
  Mcolor <- c("#8dd3c7","#ffffb3", "white")
  names(Mcolor) <- c("M0","M1",NA) #类型颜色
  Ncolor <- c("#ffffb3", "#42b5407f","#0099B47F","#00468B7F", "white")
  names(Ncolor) <- c("N0","N1","N2","N3",NA) #类型颜色
  stagecolor <- c("#bfd7ea", "#4E80AB","#a7c957", "#6a994e","white")
  names(stagecolor) <- c("stage i","stage ii","stage iii","stage iv",NA) #类型颜色
  gradecolor <- c("#bbbdc0", "#f6ddb4","#f6a34a", "white")
  names(gradecolor) <- c("G1","G2","G3",NA) #类型颜色
  MSIcolor <- c("#925fa7", "#c5b3d1", "#bbbdc0")
  names(MSIcolor) <- c("MSI-H","MSI-L","MSS") #类型颜色
  # ann_colors <- list(Cluster=clustercolor, Age= Agecolor, MSI=MSIcolor,Gender=gendercolor, Stage=stagecolor,Grade=gradecolor,
  #                    T=Tcolor,N=Ncolor,M=Mcolor) #颜色设置
  ann_colors <- list("Cluster"=clustercolor, 
                     Age= Agecolor, 
                     "Gender"=gendercolor, "Stage"=stagecolor,"Grade"=gradecolor,
                     "T"=Tcolor,"N"=Ncolor,"M"=Mcolor) #颜色设置
  colnames(clinical)=str_to_title(colnames(clinical))
  # colnames(clinical)[which(colnames(clinical)=="Msi")]="MSI"
  # annotation_col = clinical[,c("Cluster","Age","MSI","Gender","T","N","M","Stage","Grade")]
  annotation_col = clinical[,c("Cluster",
                               "Age",
                               "Gender","T","N","M","Stage","Grade")]
  
  row.names(annotation_col) <- clinical$Sample
  #在画图的时候，heatmap函数中多添加annotation_colors即可：
  p=pheatmap(as.matrix(data),cluster_rows = T,cluster_cols = T,
             color=colorRampPalette(c("navy","white","firebrick3"))(100),
             show_colnames = F,border_color = NA,scale = "row",show_rownames =T,
             annotation_col = annotation_col,# annotation_row = annotation_row,
             annotation_colors = ann_colors)
  p1=pheatmap(data,cluster_rows =T,cluster_cols = F,
              color=colorRampPalette(c("navy","white","firebrick3"))(100),
              show_colnames = F,border_color = NA,scale = "row",show_rownames =T,
              annotation_col = annotation_col,# annotation_row = annotation_row,
              annotation_colors = ann_colors)
  p_list <- list(p,p1)
  return(p_list)
}
p_list=pheatmap_clincial(train_score,TF,clinical,train_cluster,network1)
p_list_test=pheatmap_clincial(test_score,TF,clinical,test_cluster,network1)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/10.clinical")
pdf("clinical_data_pheatmap.pdf",width=8,height=6.5)
p_list[[1]]
dev.off()
pdf("clinical_data_pheatmap_no_cluster.pdf",width=8,height=6.5)
p_list[[2]]
dev.off()
pdf("test/clinical_data_pheatmap.pdf",width=8,height=6.5)
p_list_test[[1]]
dev.off()
pdf("test/clinical_data_pheatmap_no_cluster.pdf",width=8,height=6.5)
p_list_test[[2]]
dev.off()
################