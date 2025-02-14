
#####------药物治疗--------
##---准备TIDE数据库的输入数据---
#https://mp.weixin.qq.com/s/Lmu_rYr2N8hM-1U_Zkk43g （使用R编程实现参考）
#https://www.jingege.wang/2022/07/19/tide%e6%95%b0%e6%8d%ae%e6%a0%87%e5%87%86%e5%8c%96%ef%bc%9ar%e8%af%ad%e8%a8%80%e6%af%8f%e8%a1%8c%e5%87%8f%e5%8e%bb%e5%b9%b3%e5%9d%87%e5%80%bc/
#https://www.jingege.wang/2021/10/02/tide%e7%ae%97%e6%b3%95%e6%9d%a5%e9%a2%84%e6%b5%8b%e5%85%8d%e7%96%ab%e6%b2%bb%e7%96%97%e7%96%97%e6%95%88/
##---TIDE数据库输入的矩阵是标准化矩阵 这个标准化矩阵就是将基因的表达值减去该基因在正常样本的表达的均值---
rm(list=ls())
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA"
data=read.table(file.path(data_path,"TCGA-STAD-mRNA-tpm_log.txt"), header=T, row.names=1,sep="\t",check.names = F,quote="")
sample_type=read.table(file.path(data_path,"sample_type.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
data_case=data[,sample_type$sample[sample_type$type=="Tumor"]]
data_control=data[,sample_type$sample[sample_type$type=="Normal"]]
mean_control=as.data.frame(apply(data_control,1,mean))
mean_control=data.frame(gene=rownames(mean_control),value=mean_control[,1])
data_case_normalized=data_case-mean_control$value
setwd("D:/project/07.STAD_subtype/analysis0507_v7/15.drug/TIDE")
write.table(data_case_normalized,"data_case_normalized_for_TIDE.txt",sep="\t",quote=FALSE,row.names = T)
####------比较亚型的得分-----
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(cowplot)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/15.drug/TIDE")
score=read.csv("TIDE_response.csv",header=T,stringsAsFactors = F)
colnames(score)[1]="sample"
score$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", score$sample)
#score=score[,c("sample","Responder","TIDE","Dysfunction","Exclusion")]
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
train=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
train$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3-\\4", train$sample)
rownames(train)=train$sample
train$cluster=factor(train$cluster,levels=c("IM","EP"))
test$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3-\\4", test$sample)
rownames(test)=test$sample
test$cluster=factor(test$cluster,levels=c("IM","EP"))
##--连续的临床特征和亚型的关系--#
boxplot_plot <- function(data,score,character,title){
  data1=merge(data,score,by="sample")
  data1=data1[,c(character,"cluster")]
  colnames(data1)[1]="score"
  data1$cluster=factor(data1$cluster,levels=c("IM","EP"))
  cluster1_num=length(which(data1$cluster=="IM"))
  cluster2_num=length(which(data1$cluster=="EP"))
  p <- ggplot(data1,aes(cluster,score))+
    # 提琴图：
    geom_violin(aes(fill=cluster),color="white")+
    # 内部箱线图：
    geom_boxplot(fill="#a6a7ac",color="#a6a7ac",
                 width = 0.1,outlier.shape = NA)+
    scale_color_manual(values = c("#61b3de","#ffee93"))+
    scale_fill_manual(values = c("#61b3de","#ffee93"))+
    labs(x="",y=paste0(character," Score"),title=title)+
    theme_bw()+
    stat_compare_means(aes(group = cluster), 
                       method = "wilcox.test",label.x=1.25)+stat_summary(fun= mean, geom = "point",
                                                                         shape = 19, size = 2, color = "black")+
    scale_x_discrete(labels = c(paste0("IM","\n","n=(",cluster1_num,")"),paste0("EP","\n","n=(",cluster2_num,")")))+
    # geom_hline(aes(yintercept=0.6), colour="#565354", linetype="dashed")+
    # geom_hline(aes(yintercept=0.3), colour="#565354", linetype="dashed")+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.y = element_text(size = 12,
                                     face="bold"),
          legend.title=element_blank(),
          legend.position = 'none')
  
  return(p)
}
character=colnames(score)[-c(1,2,3,10)]
for(i in 1:length(character)){
  train_compare <-  boxplot_plot(train,score,character[i],"Train Datasets")
  test_compare <-  boxplot_plot(test,score,character[i],"Test Datasets")
  p_list <- list(train_compare,test_compare)
  pdf(paste0(character[i],"_train.pdf"),height=4.5,width=3)
  print(train_compare)
  dev.off()
  pdf(paste0(character[i],"_test.pdf"),height=4.5,width=3)
  print(test_compare)
  dev.off()
  pdf(paste0(character[i],"_train_test.pdf"),height=4.5,width=6)
  print(plot_grid(plotlist=p_list,ncol=2))
  dev.off()
}

##---TIDE疗效和cluster的关系------
clinical_barplot <- function(cluster,score,character,title){
  data=merge(cluster,score,by="sample")
  data1=data[,c("cluster",character)]
  colnames(data1)[2]="clinical"
  table_clinical=table(data1[,"clinical"],data1[,'cluster'])
  fisher_res=fisher.test(table_clinical)
  pv=round(fisher_res$p.value,3)
  if(pv<0.001){pv="P<0.001"}else{pv=paste("P=",round(pv,3),sep="")}
  color=c("#ffee93", "#bbbdc0")
  data1$clinical=factor(data1$clinical,levels=c("Response","None-response"))
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
score$Responder[score$Responder=="True"]="Response"
score$Responder[score$Responder=="False"]="None-response"
p1_train=clinical_barplot(train,score,"Responder","Train Datasets")
p2_test=clinical_barplot(test,score,"Responder","Test Datasets")
p_list=list(p1_train,p2_test)
pdf("response_train_test.pdf",height=4.5,width=6)
print(plot_grid(plotlist=p_list,ncol=2))
dev.off()
pdf("response_train.pdf",height=4.5,width=4)
print(p1_train)
dev.off()
pdf("response_test.pdf",height=4.5,width=4)
print(p2_test)
dev.off()
##---化疗药物--------
#---不同亚型的药物预测 oncopredict-----
##参考https://www.jingege.wang/2022/09/21/%e8%8d%af%e7%89%a9%e6%95%8f%e6%84%9f%e6%80%a7%e9%a2%84%e6%b5%8br%e5%8c%85%e4%b9%8boncopredict/
##--oncopredict 得分越高敏感性越好 越适合--#
#install.packages("oncoPredict")
rm(list=ls())
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/15.drug/oncoPredict")
GDSC2_Expr = readRDS('GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_Res = readRDS("GDSC2_Res.rds")
GDSC2_Res <- exp(GDSC2_Res)
######读入你需要做预测的表达量矩阵
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA"
data=read.table(file.path(data_path,"TCGA-STAD-mRNA-tpm_log_case.txt"), header=T, row.names=1,sep="\t",check.names = F,quote="")
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(data))
setwd("D:/project/07.STAD_subtype/analysis0507_v7/15.drug/oncoPredict")
dataset=as.matrix(data)
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = dataset,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )
####------比较亚型的oncopredict得分-----
####------比较亚型的oncopredict得分-----
##--oncopredict 得分越低约好--#
##---选择IC50得分低于10的药物---
rm(list=ls())
library(stringr)
library(ComplexHeatmap)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/15.drug/oncoPredict")
score=read.csv("calcPhenotype_Output/DrugPredictions.csv",header=T,stringsAsFactors = F,check.names = F)
colnames(score)[1]="sample"
colnames(score)=str_split_fixed(colnames(score),"_",2)[,1]
rownames(score)=score$sample
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
train=read.table(file.path(cluster_path,"train_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table(file.path(cluster_path,"test_cluster_rename.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
train$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3-\\4", train$sample)
rownames(train)=train$sample
train$cluster=factor(train$cluster,levels=c("IM","EP"))
test$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3-\\4", test$sample)
rownames(test)=test$sample
test$cluster=factor(test$cluster,levels=c("IM","EP"))
data=merge(train,score,by="sample")
rownames(data)=data$sample
data1=data[,-c(1,2)]
data_case=as.data.frame(t(data1[data$cluster=="IM",]))
data_control=as.data.frame(t(data1[data$cluster=="EP",]))
case_max=apply(data_case,1,max)
case_mean=apply(data_case,1,mean)
control_max=apply(data_control,1,max)
control_mean=apply(data_control,1,mean)
fc=case_mean/control_mean
de=round(case_mean-control_mean,2)
fc_de=data.frame(drug=rownames(data_case),fc=fc,de=de)
sigDrug=NULL
result=NULL
for(i in 3:dim(data)[2]){
  # if(sd(data[,i])<0.001){next}
  fit=wilcox.test(data[,i] ~ data[,2])
  pvalue=fit$p.value
  temp=data.frame(drug=colnames(data)[i],p=round(pvalue,3))
  result=rbind(result,temp)
  if(pvalue<0.05){
    sigDrug=c(sigDrug, colnames(data)[i])
  }
}
fc_de=fc_de[order(abs(fc_de$de),fc_de$fc,decreasing = T),]
result=merge(result,fc_de,by="drug")
sig=fc_de[fc_de$drug %in% sigDrug,]
sig= sig[sig$drug %in% Reduce(intersect,list(rownames(data_case)[which(case_mean<10|control_mean<10)],rownames(data_case)[which(case_max<10|control_max<10)],sig$drug)),]
# sig= sig[sig$drug %in% intersect(rownames(data_case)[which(case_mean<10|control_mean<10)],sig$drug),]
write.table(result,"drug_p_de.xls",sep="\t",quote=FALSE,row.names = F)
write.table(sig,"sig_drug.xls",sep="\t",quote=FALSE,row.names = F)
cluster=data[,c(1,2)]
cluster$cluster=factor(cluster$cluster,levels=c("IM","EP"))
cluster=cluster[order(cluster$cluster),]
data2=cbind(data_case,data_control)
data2=merge(data2,sig,by.x="row.names",by.y="drug")
rownames(data2)=data2[,1]
data2=data2[,-1]
data2=data2[-c(7,10),]
Up=data2%>%dplyr::filter(de>0)
Down=data2%>%dplyr::filter(de<0)
normal=data2%>%dplyr::filter(de==0)
expr=rbind(Up,Down,normal)
expr=expr[,cluster$sample]
drug=c(rep("EP",nrow(Up)),rep("IM",nrow(Down)),rep("both",nrow(normal)))
scaled_expr=t(scale(t(expr))) 
ha1=HeatmapAnnotation(df=data.frame(
  cluster=cluster$cluster),
  col=list(cluster=c("IM"= "#61b3de", "EP"="#ffee93")))
ha2=rowAnnotation(df=data.frame(
  drug=c(rep("EP",nrow(Up)),rep("IM",nrow(Down)),rep("both",nrow(normal)))),
  col=list(drug=c("IM"="#61b3de","EP"="#ffee93","both"="#E18727FF")))
library(circlize)
names=str_split_fixed(rownames(scaled_expr),"[.]",2)[,1]
pdf(file="complexheatmap_nocluster.pdf",width=6.5,height=7)
Heatmap(scaled_expr,name="Estimated IC50",top_annotation=ha1,left_annotation =ha2 ,
        cluster_rows=TRUE,cluster_columns=FALSE,show_row_names=TRUE,
        show_column_names=FALSE,
        column_split = cluster$cluster,
        row_split = drug,
        #row_labels = names,
        #column_names_gp=gpar(fontsize=5),
        row_names_gp = gpar(fontsize=8))
dev.off()

##--显著的药物箱式图--#
library(ggplot2)
library(ggpubr)
library(cowplot)
boxplot_plot <- function(data1,character){
  data1=data1[,c(character,"cluster")]
  colnames(data1)[1]="score"
  data1$cluster=factor(data1$cluster,levels=c("IM","EP"))
  cluster1_num=length(which(data1$cluster=="IM"))
  cluster2_num=length(which(data1$cluster=="EP"))
  drug=str_split_fixed(character,"[.]",2)[1]
  p <- ggplot(data1,aes(x=cluster,y=score,fill=cluster))+geom_boxplot(width=0.8,position=position_dodge(0.9),
                                                                      #outlier.colour = NA,
                                                                      notch = T)+
    labs(y=paste0("Estimated IC50 of ",character),x="")+scale_fill_manual(values = c("#61b3de", "#ffee93"))+
    stat_compare_means(aes(group = cluster), 
                       method = "wilcox.test",label.x=1)+
    theme_bw()+
    #theme(axis.text.x = element_text(angle = 0,vjust = 0.85,hjust = 0.75))+
    theme(legend.position = 'none')+
    scale_x_discrete(labels = c(paste0("IM","\n","n=(",cluster1_num,")"),paste0("EP","\n","n=(",cluster2_num,")")))+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.y = element_text(size = 12,
                                     face="bold"))
  pdf(paste0(drug,"_divide.pdf"),width=3,height=5)
  print(p)
  dev.off()
  return(p)
}
drug=rownames(scaled_expr)
p_list=list()
for(i in 1:length(drug)){
  p_list[[i]]=boxplot_plot(data,drug[i])
}
pdf("sig_drug.pdf",width=9,height=18)
plot_grid(plotlist = p_list,ncol=4)
dev.off()
###############