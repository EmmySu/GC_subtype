##------亚型的功能---------#
##----得到两个亚群差异的甲基化位点-------
rm(list=ls())
library(reshape2)
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
cluster=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
rownames(cluster)=cluster$sample
cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
meth_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/meth"
data=readRDS(file.path(meth_path,"train_data.rds"))
IM=data[,cluster$sample[cluster$cluster=="IM"]]
EP=data[,cluster$sample[cluster$cluster=="EP"]]
##--获得差异甲基化位点---#
get_diff_point <- function(data1,data2){
  mean_data1=apply(data1,1,mean)
  mean_data2=apply(data2,1,mean)
  chazhi=mean_data2-mean_data1
  #sig_chazhi=rownames(data1)[which(abs(chazhi)>=0.2)]
  data=cbind(data1,data2)
  data=as.data.frame(data)
  data$chazhi=chazhi
  data$p=0
  for(i in 1:dim(data)[1]){
    p=t.test(data2[i,],data1[i,])$p.value
    data$p[i]=p
  }
  data$padj=p.adjust(data$p,method="fdr",length(data$p))
  return(data)
}
IMvsEP=get_diff_point(EP,IM)
sig_T_N=rownames(IMvsEP)[which(abs(IMvsEP$chazhi)>=0.2 & IMvsEP$padj<0.05)]#
sig_T_N2=rownames(IMvsEP)[which(abs(IMvsEP$chazhi)>=0.1 & IMvsEP$padj<0.05)]#
sig=IMvsEP[sig_T_N2,]###
fc_chaizhi=IMvsEP[,c((dim(IMvsEP)[2]-2) : dim(IMvsEP)[2])]
fc_chaizhi=cbind(ID=rownames(fc_chaizhi),fc_chaizhi)
diff_meth=list(IMvsEP=IMvsEP,sig_T_N=sig_T_N2,sig=sig,fc_chaizhi=fc_chaizhi)
dir.create("D:/project/07.STAD_subtype/analysis0507_v7/14.function/meth")
setwd("D:/project/07.STAD_subtype/analysis0507_v7/14.function/meth")
saveRDS(IMvsEP,"IMvsEP.rds")
saveRDS(diff_meth,"diff_meth_result.rds")
write.table(sig,"sig_meth.xls",sep="\t",quote=FALSE,row.names=T)
write.table(fc_chaizhi,"fc_chazhi.xls",sep="\t",quote=FALSE,row.names=F)
##---看差异的甲基化位点都是哪些基因---#
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/14.function/meth")
sig=read.table("sig_meth.xls",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
platform_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/meth"
id=read.table(file.path(platform_path,"promoter_cg_hg38.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
sig_gene=id[id$cg %in% rownames(sig),]
gene=unique(sig_gene$gene)##2311
aa=as.data.frame(table(sig_gene$gene))
write.table(gene,"sig_cg_gene.xls",sep="\t",quote=FALSE,row.names=F)
write.table(sig_gene,"sig_cg_gene_point.xls",sep="\t",quote=FALSE,row.names=F)
fc_chaizhi=read.table("fc_chazhi.xls",sep="\t",header=T,stringsAsFactors = F,quote="")
gene_set_path="D:/project/07.STAD_subtype/analysis0507_v7/05.cox_TF_cox_CpG_cor"
##读取背景基因集合
gene_set=read.table(file.path(gene_set_path,"TF_CpG_cor_sig.xls"),header=T,sep="\t",stringsAsFactors = F)
cg=unique(gene_set$cg)
sig_cg=intersect(rownames(sig),cg)
p_cg=fc_chaizhi[fc_chaizhi$ID %in% cg,]
length(which(p_cg$padj<0.05))
##----差异甲基化位点热图--------
rm(list=ls())
library(magrittr)
library(ComplexHeatmap)
data_path="D:/project/07.STAD_subtype/analysis0507_v7//00.data/meth"
setwd("D:/project/07.STAD_subtype/analysis0507_v7/14.function/meth")
data=readRDS("diff_meth_result.rds")
IMvsEP=data$IMvsEP
sig=data$sig
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
cluster=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
rownames(cluster)=cluster$sample
cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
cluster=cluster[order(cluster$cluster),]
Up=sig%>%dplyr::filter(chazhi >=0.1)
Down=sig%>%dplyr::filter(chazhi <= -0.1)
expr=rbind(Up,Down)
expr=expr[,cluster$sample]
###---标准化----###
scaled_expr=t(scale(t(expr),
                    center = T, scale=T)) 
ha1=HeatmapAnnotation(df=data.frame(
  cluster=cluster$cluster),
  col=list(cluster=c("IM"= "#61b3de", "EP"="#ffee93")))
ha2=rowAnnotation(df=data.frame(
  CpG=c(rep("Up",nrow(Up)),rep("Down",nrow(Down)))),
  col=list(CpG=c("Up"="tomato","Down"="steelblue")))
library(circlize)
pdf(file="complexheatmap_nocluster.pdf",width=5,height=5)
Heatmap(scaled_expr,name="CpG",top_annotation=ha1,left_annotation =ha2 ,
        cluster_rows=FALSE,cluster_columns=FALSE,show_row_names=FALSE,
        show_column_names=FALSE,
        column_names_gp=gpar(fontsize=12))
dev.off()
pdf(file="complexheatmap_cluster.pdf",width=5,height=5)
Heatmap(scaled_expr,name="CpG",top_annotation=ha1,left_annotation =ha2 ,
        cluster_rows=FALSE,cluster_columns=TRUE,show_row_names=FALSE,
        show_column_names=FALSE,
        column_names_gp=gpar(fontsize=12))
dev.off()
library(pheatmap)
annotation_col = data.frame( cluster = cluster$cluster)
rownames(annotation_col) = cluster$sample
annotation_row = data.frame( CpG=c(rep("Up",nrow(Up)),rep("Down",nrow(Down))))
rownames(annotation_row) = c(rownames(Up),rownames(Down))
ann_colors = list( cluster = c(IM ="#61b3de", EP = "#ffee93"),
                   CpG=c(Up="#61b3de",Down="#ffee93"))
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
         show_rownames = F,
         color = colorRampPalette(rev(c("#BC3C29FF","#FCFCE3","#0072B5FF")))(length(bk)),
         fontsize_number =8, number_color = "white",
         gaps_col =cluster1.length,
         gaps_row = up_length,
         annotation_row = annotation_row,
         
         annotation_col =annotation_col, 
         annotation_colors =ann_colors)
dev.off()

##----得到两个亚群差异的mRNA-------
rm(list=ls())
library(reshape2)
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
cluster=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", cluster$sample)
rownames(cluster)=cluster$sample
cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA"
data=read.table(file.path(data_path,"TCGA-STAD-mRNA-tpm_log_case.txt"), header=T, row.names=1,sep="\t",check.names = F,quote="")
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
sample=intersect(cluster$sample,colnames(data))
cluster=cluster[sample,]
cluster=cluster[order(cluster$cluster),]
IM=data[,cluster$sample[cluster$cluster=="IM"]]
EP=data[,cluster$sample[cluster$cluster=="EP"]]
##--获得差异mRNA---#
get_diff_point <- function(data1,data2){
  mean_data1=apply(data1,1,mean)
  mean_data2=apply(data2,1,mean)
  log2FC=mean_data2-mean_data1
  #sig_log2FC=rownames(data1)[which(abs(log2FC)>=2)]
  data=cbind(data1,data2)
  data=as.data.frame(data)
  data$log2FC=log2FC
  data$p=0
  for(i in 1:dim(data)[1]){
    p=t.test(data2[i,],data1[i,])$p.value
    data$p[i]=p
  }
  data$padj=p.adjust(data$p,method="fdr",length(data$p))
  return(data)
}
IMvsEP=get_diff_point(EP,IM)
sig_T_N=rownames(IMvsEP)[which(abs(IMvsEP$log2FC)>=1 & IMvsEP$padj<0.05)]#
sig_T_N2=rownames(IMvsEP)[which(abs(IMvsEP$log2FC)>=log2(1.5) & IMvsEP$padj<0.05)]#
sig=IMvsEP[sig_T_N2,]###
fc_p=IMvsEP[,c((dim(IMvsEP)[2]-2) : dim(IMvsEP)[2])]
fc_p=cbind(ID=rownames(fc_p),fc_p)
diff_mRNA=list(IMvsEP=IMvsEP,sig_T_N=sig_T_N2,sig=sig,fc_p=fc_p)
dir.create("D:/project/07.STAD_subtype/analysis0507_v7/14.function/mRNA")
setwd("D:/project/07.STAD_subtype/analysis0507_v7/14.function/mRNA")
saveRDS(IMvsEP,"IMvsEP.rds")
saveRDS(diff_mRNA,"diff_mRNA_result.rds")
write.table(sig,"sig_mRNA.xls",sep="\t",quote=FALSE,row.names=T)
write.table(fc_p,"log2FC_p.xls",sep="\t",quote=FALSE,row.names=F)
gene_set_path="D:/project/07.STAD_subtype/analysis0507_v7/05.cox_TF_cox_CpG_cor"
gene_set=read.table(file.path(gene_set_path,"TF_CpG_cor_sig.xls"),header=T,sep="\t",stringsAsFactors = F)
TF=unique(gene_set$TF)
intersect(TF,rownames(sig))
##---绘制差异的mRNA的热图----
rm(list=ls())
library(magrittr)
library(ComplexHeatmap)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/14.function/mRNA")
data=read.delim("sig_mRNA.xls",sep="\t",header=T,stringsAsFactors=F,check.names=F)
data=data[order(abs(data$log2FC)),]
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
cluster=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", cluster$sample)
rownames(cluster)=cluster$sample
cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
sample=intersect(cluster$sample,colnames(data))
cluster=cluster[sample,]
cluster=cluster[order(cluster$cluster),]
Up=data%>%dplyr::filter(log2FC >=log2(1.5))
Down=data%>%dplyr::filter(log2FC <= log2(1/1.5))
expr=rbind(Up,Down)
expr=expr[,cluster$sample]
###---标准化----###
scaled_expr=t(scale(t(expr),
                    center = T, scale=T)) 
ha1=HeatmapAnnotation(df=data.frame(
  cluster=cluster$cluster),
  col=list(cluster=c("IM"= "#61b3de", "EP"="#ffee93")))
ha2=rowAnnotation(df=data.frame(
  mRNA=c(rep("Up",nrow(Up)),rep("Down",nrow(Down)))),
  col=list(mRNA=c("Up"="#61b3de","Down"="#ffee93")))
library(circlize)
pdf(file="complexheatmap_nocluster.pdf",width=5,height=5)
Heatmap(scaled_expr,name="mRNA",top_annotation=ha1,left_annotation =ha2 ,
        cluster_rows=FALSE,cluster_columns=FALSE,show_row_names=FALSE,
        show_column_names=FALSE,
        column_names_gp=gpar(fontsize=12))
dev.off()
pdf(file="complexheatmap_cluster.pdf",width=5,height=5)
Heatmap(scaled_expr,name="mRNA",top_annotation=ha1,left_annotation =ha2 ,
        cluster_rows=FALSE,cluster_columns=TRUE,show_row_names=FALSE,
        show_column_names=FALSE,
        column_names_gp=gpar(fontsize=12))
dev.off()
library(pheatmap)
annotation_col = data.frame( cluster = cluster$cluster)
rownames(annotation_col) = cluster$sample
annotation_row = data.frame( mRNA=c(rep("Up",nrow(Up)),rep("Down",nrow(Down))))
rownames(annotation_row) = c(rownames(Up),rownames(Down))
ann_colors = list( cluster = c(IM ="#61b3de", EP = "#ffee93"),
                   mRNA=c(Up="#61b3de",Down="#ffee93"))
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
         show_rownames = F,
         color = colorRampPalette(rev(c("#BC3C29FF","#FCFCE3","#0072B5FF")))(100),
         fontsize_number =8, number_color = "white",
         gaps_col =cluster1.length,
         gaps_row = up_length,
         annotation_row = annotation_row,
         
         annotation_col =annotation_col, 
         annotation_colors =ann_colors)
dev.off()
#---cluster difff huoshantu 和标记显著的TF-----
rm(list=ls())
library(ggplot2)
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA"
setwd("D:/project/07.STAD_subtype/analysis0507_v7/14.function/mRNA")
data=readRDS("IMvsEP.rds")
data=data[,c("log2FC","p","padj")]
colnames(data)[1]=c("log2FoldChange")
sig=read.table("sig_mRNA.xls",sep="\t",header=T,stringsAsFactors = F,quote="")
TF_path="D:/project/07.STAD_subtype/analysis0507_v7/04.OS/TF"
TF=read.table(file.path(TF_path,"single_cox_sig.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
TF_sig=unique(TF$Characteristics)
data$label=NA
TF_sig=intersect(rownames(sig),TF_sig)
data[TF_sig,4]=TF_sig
library(ggplot2)
library(ggrepel)
pdf("cluster_diff_huoshantu.pdf",width=5,height=5)
ggplot(data,aes(log2FoldChange, -log10(padj)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(log2(1/1.5),log2(1.5)), linetype = "dashed", color = "#999999")+
  geom_point(aes(size=-log10(padj), color= -log10(padj)))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  scale_size_continuous(range = c(1,2))+
  theme_bw()+
  theme(panel.grid = element_blank()) +
  geom_text_repel(data = data, aes(x =log2FoldChange, 
                                   y =-log10(padj), 
                                   label =label),
                  size = 2,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30) ,
                  show.legend = FALSE)
dev.off()
########
##---两个亚型差异的mRNA做kegg功能富集分析------
rm(list=ls())
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(ggplot2)
library(cols4all)
##---function------
##---ID转化--#
ID_trans <- function(data){
  symbol=rownames(data)
  keyword <- as.character(symbol)
  res <-mapIds(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
               keys=keyword,
               column="ENTREZID" , #clolumns参数是你要转换的ID类型是什么，只能选择一个。
               keytype="SYMBOL")#函数里面的keytype与keys参数是对应的，keys是你输入的那些数据，keytype是指这些数据是属于什么类型的数据。 
  res <- as.data.frame(res)
  res$gene <- rownames(res);colnames(res)[1] <- "ENTREZID"
  temp <- res[!is.na(res$ENTREZID),]
  temp <- temp[,c(2,1)]
  colnames(temp) <- c("symbol","gene.id")
  return(temp)
}
##--KEGG通路富集--#
KEGG_Enrichment <- function(temp){
  gene=as.character(temp$gene.id)
  ekk <- enrichKEGG(gene=gene,organism = "hsa",pAdjustMethod = "BH",pvalueCutoff=0.05)
  ekk <- as.data.frame(ekk)
  if(dim(ekk)[1]>0){
    newdata <- ekk[order(ekk$Count, ekk$p.adjust),]
    newdata$Description = factor(newdata$Description,levels = newdata$Description)
    newdata$Pop_Hits <- sapply(newdata$BgRatio, function(x) {nn <- unlist(strsplit(x,split="[/]"));mm <- nn[1];return(mm)})
    newdata$Pop_Hits <- as.numeric(newdata$Pop_Hits)
  }else{
    newdata=NULL
  }
  
  return(newdata)
}
##--barplot---##
enrich_barplot <- function(dt,title){
  dt=dt[,c("Description","-log10pvalue")]
  dt$group <- case_when(dt$`-log10pvalue` > 0 ~ 'up',
                        dt$`-log10pvalue` < 0 ~ 'down')
  #指定因子；调整顺序：
  dt$Description=factor(dt$Description,levels=rev(dt$Description))
  p <- ggplot(dt,
              aes(x =`-log10pvalue`, y = Description, fill = `-log10pvalue`)) + #数据映射
    geom_col(width = 0.5) + #绘制添加条形图
    #scale_x_discrete(labels=function(x) str_wrap(x, width=50))+###太长的通路名字会换行
    theme_bw()
  #自定义主题调整：
  mytheme <- theme(
   # legend.position = 'none',
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = 'grey60',size = 1.1),
    axis.text = element_text(size = 12)
  )
  p1 <- p + mytheme
  #先根据上下调标签拆分数据框：
  up <- dt[which(dt$`-log10pvalue` > 0),]
  down <- dt[which(dt$`-log10pvalue` < 0),]
  #添加上调pathway标签：
  p2 <- p1 +
    geom_text(data = up,
              aes(x = -0.2, y = Description, label = Description),
              size = 3.5,
              hjust = 1) #标签右对齐
  #添加下调pathway标签：
  p3 <- p2 +
    geom_text(data = down,
              aes(x = 0.2, y = Description, label = Description),
              size = 3.5,
              hjust = 0) #标签左对齐
  p4 <- p3 +
    #scale_x_continuous(breaks=seq(-4, 6, 2)) + #x轴刻度修改
    labs(x = '-log(Pvalue)', y = ' ', title = title) + #修改x/y轴标签、标题添加
    theme(plot.title = element_text(hjust = 0.5, size = 14)) #主标题居中、字号调整
  mycol <- c("#79C178","#356085")
  p5 <- p4 +
    #scale_fill_manual(values = mycol)
    scale_color_continuous_c4a_div('sunset', mid = 0, reverse = F) +
    scale_fill_continuous_c4a_div('sunset', mid = 0, reverse = F) 
}
enrich_lolliplot <- function(dt,title){
  dt=dt[,c("Description","Count","-log10pvalue")]
  dt$group <- case_when(dt$`-log10pvalue` > 0 ~ 'up',
                        dt$`-log10pvalue` < 0 ~ 'down')
  #指定因子；调整顺序：
  dt$Description=factor(dt$Description,levels=rev(dt$Description))
  p <- ggplot(dt,
              aes(x =`-log10pvalue`, y = Description, fill = `-log10pvalue`)) + #数据映射
    geom_col(aes(fill = `-log10pvalue`), width = 0.1) + #绘制添加条形图
    geom_point(aes(size = Count,
                   color = `-log10pvalue`)) +
    scale_size_continuous(range = c(1, 3)) +
    #scale_x_discrete(labels=function(x) str_wrap(x, width=50))+###太长的通路名字会换行
    theme_bw()
  #自定义主题调整：
  mytheme <- theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = 'grey60',size = 1.1),
    axis.text = element_text(size = 10)
  )
  p1 <- p + mytheme
  #先根据上下调标签拆分数据框：
  up <- dt[which(dt$`-log10pvalue` > 0),]
  down <- dt[which(dt$`-log10pvalue` < 0),]
  #添加上调pathway标签：
  p2 <- p1 +
    geom_text(data = up,
              aes(x = -0.2, y = Description, label = Description),
              size = 3.5,
              hjust = 1) #标签右对齐
  #添加下调pathway标签：
  p3 <- p2 +
    geom_text(data = down,
              aes(x = 0.2, y = Description, label = Description),
              size = 3.5,
              hjust = 0) #标签左对齐
  p4 <- p3 +
    #scale_x_continuous(breaks=seq(-4, 6, 2)) + #x轴刻度修改
    labs(x = '-log(Pvalue)', y = ' ', title = title) + #修改x/y轴标签、标题添加
    theme(plot.title = element_text(hjust = 0.5, size = 14)) #主标题居中、字号调整
  mycol <- c("#79C178","#356085")
  p5 <- p4 +
    #scale_fill_manual(values = mycol)
    scale_color_continuous_c4a_div('sunset', mid = 0, reverse = F) +
    scale_fill_continuous_c4a_div('sunset', mid = 0, reverse = F) + scale_x_continuous(breaks = seq(-20,20,by = 10),
                                                                                       labels = abs(seq(-20,20,by = 10))) 
}
setwd("D:/project/07.STAD_subtype/analysis0507_v7/14.function/mRNA")
data=read.delim("sig_mRNA.xls",sep="\t",header=T,stringsAsFactors=F,check.names=F)
data=data[order(abs(data$log2FC)),]
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
cluster=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
cluster$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", cluster$sample)
rownames(cluster)=cluster$sample
cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
sample=intersect(cluster$sample,colnames(data))
cluster=cluster[sample,]
cluster=cluster[order(cluster$cluster),]
Up=data%>%dplyr::filter(log2FC >=log2(1.5))
Down=data%>%dplyr::filter(log2FC <= log2(1/1.5))


up_gene_id <- ID_trans(Up)
down_gene_id <- ID_trans(Down)
up_kegg <- KEGG_Enrichment(up_gene_id)
down_kegg <- KEGG_Enrichment(down_gene_id)
up_kegg$`-log10pvalue`=-log10(up_kegg$pvalue)
up_kegg=up_kegg[order(up_kegg$`-log10pvalue`,decreasing = T),]
down_kegg$`-log10pvalue`=log10(down_kegg$pvalue)
down_kegg=down_kegg[order(down_kegg$`-log10pvalue`,decreasing = T),]
dt=rbind(up_kegg,down_kegg)
write.table(dt,"kegg_enrich.xls",sep="\t",quote=FALSE,row.names = F)
p=enrich_barplot(dt,"IM vs EP Enriched KEGG Pathway")
ggsave("kegg_barplot.pdf",p, height =7.5, width = 12.5)
p=enrich_lolliplot(dt,"IM vs EP Enriched KEGG Pathway")
ggsave("kegg_lolliplot.pdf",p, height =8.5, width = 10.5)
