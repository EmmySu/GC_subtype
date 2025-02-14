#------一.TCGA STAD mRNA数据提取、log 转化和获得样本类型----------
##---下载数据----
library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
library(reshape2)
projects <- getGDCprojects()
projects <- projects %>% 
  as.data.frame() %>% 
  select(project_id,tumor) %>% 
  filter(grepl(pattern="TCGA",project_id))
###-----function----#
get_expression <- function(cancer){
  ## 1.查询信息
  query.exp=GDCquery(project=cancer, 
                     data.category="Transcriptome Profiling",
                     data.type="Gene Expression Quantification",
                     workflow.type="STAR - Counts")
  ## 2.正式下载
  GDCdownload(query.exp)
  ## 3.多个数据合并
  pre.exp=GDCprepare(query=query.exp)
  ## 4.提取表达量数据
  countsdata=SummarizedExperiment::assay(pre.exp,1)
  fpkmdata=SummarizedExperiment::assay(pre.exp,5)
  tpmdata=SummarizedExperiment::assay(pre.exp,4)
  gene_id=data.frame(id=rowData(pre.exp)@listData[["gene_id"]], gene_name= rowData(pre.exp)@listData[["gene_name"]],gene_type=rowData(pre.exp)@listData[["gene_type"]])
  counts=cbind(gene_id,countsdata)
  fpkm=cbind(gene_id,fpkmdata)
  tpm=cbind(gene_id,tpmdata)
  #临床信息
  clinical=GDCquery_clinic(project=cancer, type="clinical")
  ## 5.保存数据
  filename1=paste0(cancer,"-counts.txt")
  filename2=paste0(cancer,"-fpkm.txt")
  filename3=paste0(cancer,"-tpm.txt")
  filename4=paste0(cancer,"-clinical.txt")
  filename5=paste0(cancer,"-gene_id.txt")
  write.table(counts,filename1,sep="\t",col.names=T,row.names=F,quote=F) 
  write.table(fpkm,filename2,sep="\t",col.names=T,row.names=F,quote=F) 
  write.table(tpm,filename3,sep="\t",col.names=T,row.names=F,quote=F) 
  write.table(clinical,filename4,sep="\t",col.names=T,row.names=F,quote=F) 
  write.table(gene_id,filename5,sep="\t",col.names=T,row.names=F,quote=F) 
}
get_expression("TCGA-STAD")
##---提取mRNA count数据-------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA")
TCGA=read.table("TCGA-STAD-counts.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
bioclass=read.table("biotype_classs.txt",sep="\t",header=F,stringsAsFactors = F,quote="")
TCGA=TCGA[TCGA$gene_type %in% bioclass$V1[bioclass$V3=="gene"],]
length(unique(TCGA$gene_name))
factors=factor(TCGA$gene_name) #
temps=tapply(TCGA[,4], factors, mean)
end=dim(TCGA)[2]
for(j in 5:end){
  temp <- tapply(TCGA[,j], factors, mean)#多个探针匹配到一个基因上
  temps <- cbind(temps, temp)
}
temps=as.data.frame(temps)
colnames(temps)=colnames(TCGA)[-c(1:3)]
TCGA=cbind(GeneID=rownames(temps),temps)
write.table(TCGA,"TCGA-STAD-mRNA-counts.txt",sep="\t",quote=FALSE,row.names = F)
##---提取mRNA fpkm数据------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA")
TCGA=read.table("TCGA-STAD-fpkm.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
bioclass=read.table("biotype_classs.txt",sep="\t",header=F,stringsAsFactors = F,quote="")
TCGA=TCGA[TCGA$gene_type %in% bioclass$V1[bioclass$V3=="gene"],]
length(unique(TCGA$gene_name))
factors=factor(TCGA$gene_name) #
temps=tapply(TCGA[,4], factors, mean)
end=dim(TCGA)[2]
for(j in 5:end){
  temp <- tapply(TCGA[,j], factors, mean)#多个探针匹配到一个基因上
  temps <- cbind(temps, temp)
}
temps=as.data.frame(temps)
colnames(temps)=colnames(TCGA)[-c(1:3)]

TCGA=cbind(GeneID=rownames(temps),temps)
write.table(TCGA,"TCGA-STAD-mRNA-fpkm.txt",sep="\t",quote=FALSE,row.names = F)
###----tpm-------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA")
TCGA=read.table("TCGA-STAD-tpm.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
bioclass=read.table("biotype_classs.txt",sep="\t",header=F,stringsAsFactors = F,quote="")
TCGA=TCGA[TCGA$gene_type %in% bioclass$V1[bioclass$V3=="gene"],]
factors=factor(TCGA$gene_name) #
temps=tapply(TCGA[,4], factors, mean)
end=dim(TCGA)[2]
for(j in 5:end){
  temp <- tapply(TCGA[,j], factors, mean)#多个探针匹配到一个基因上
  temps <- cbind(temps, temp)
}
temps=as.data.frame(temps)
colnames(temps)=colnames(TCGA)[-c(1:3)]
TCGA=cbind(GeneID=rownames(temps),temps)
write.table(TCGA,"TCGA-STAD-mRNA-tpm.txt",sep="\t",quote=FALSE,row.names = F)
###---log------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA")
data=read.table("TCGA-STAD-mRNA-counts.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F,row.names = 1)
data=apply(data,2,function(x){m=log2(x+1);return(m)})
data=as.matrix(data)
write.table(data,"TCGA-STAD-mRNA-counts_log.txt",sep="\t",quote=FALSE,row.names = T)

data=read.table("TCGA-STAD-mRNA-fpkm.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F,row.names = 1)
data=apply(data,2,function(x){m=log2(x+1);return(m)})
data=as.matrix(data)
write.table(data,"TCGA-STAD-mRNA-fpkm_log.txt",sep="\t",quote=FALSE,row.names = T)

data=read.table("TCGA-STAD-mRNA-tpm.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F,row.names = 1)
data=apply(data,2,function(x){m=log2(x+1);return(m)})
data=as.matrix(data)
write.table(data,"TCGA-STAD-mRNA-tpm_log.txt",sep="\t",quote=FALSE,row.names = T)
##################
##---mRNA 样本类型-------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA")
TCGA=read.table("TCGA-STAD-mRNA-counts.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F,row.names = 1)
group=sapply(strsplit(colnames(TCGA),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
samplesNT=colnames(TCGA)[which(group==1)]
samplesTP=colnames(TCGA)[which(group==0)]
type=data.frame(sample=c(samplesNT,samplesTP),type=c(rep("Normal",length(samplesNT)),
                                                     rep("Tumor",length(samplesTP))))
write.table(type,"sample_type.txt",sep="\t",quote=FALSE,row.names = F)
##--mRNA tpm log case-----
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA")
data=read.table("TCGA-STAD-mRNA-tpm_log.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F,row.names = 1)
sample_type=read.table("sample_type.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
data1=data[,sample_type$sample[sample_type$type=="Tumor"]]
write.table(data1,"TCGA-STAD-mRNA-tpm_log_case.txt",sep="\t",quote=FALSE,row.names = T)
#################
#----二、获得TF-------
##-RcisTarget----#
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/TF")
library(RcisTarget)
data(motifAnnotations_hgnc)
data(motifAnnotations_hgnc_v9)
TF1=unique(motifAnnotations$TF)
TF2=unique(motifAnnotations_hgnc_v9$TF)
setdiff(TF2,TF1)
setdiff(TF1,TF2)
result=data.frame(TF=TF1,source="RcisTarget")
write.table(result,"TF_RcisTarget.xls",sep="\t",quote=FALSE,row.names = F)
##---union----
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/TF")
TFDB=read.table("AnimalATFDB3.0.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
TFDB=unique(TFDB)
colnames(TFDB)=c("TF","source")
hTFtarget=read.table("hTFtarget.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
hTFtarget=unique(hTFtarget)
colnames(hTFtarget)=c("TF","source")
Wei=read.table("Wei.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
Wei$source="Wei"
colnames(Wei)=c("TF","source")
Wei=unique(Wei)
TRRUST=read.table("TRRUST.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
TRRUST$source="TRRUST"
TRRUST=TRRUST[,-2]
colnames(TRRUST)=c("TF","source")
TRRUST=unique(TRRUST)
RcisTarget=read.table("TF_RcisTarget.xls",sep="\t",header=T,stringsAsFactors = F,quote="")
colnames(RcisTarget)=c("TF","source")
RcisTarget=unique(RcisTarget)
hgnc=read.table("hs_hgnc_tfs.txt",sep="\t",header=F,stringsAsFactors = F,quote="")
hgnc$source="hgnc"
colnames(hgnc)=c("TF","source")
hgnc=unique(hgnc)
TF=Reduce(union,list(hTFtarget$TF,Wei$TF,TRRUST$TF,RcisTarget$TF,TFDB$TF,hgnc$TF))
TF_inter=Reduce(intersect,list(hTFtarget$TF,Wei$TF,TRRUST$TF,RcisTarget$TF,TFDB$TF,hgnc$TF))
data=rbind(hTFtarget,Wei,TRRUST,RcisTarget,TFDB,hgnc)
data=unique(data)
count=as.data.frame(table(data$TF))
data1=data[data$TF %in% count$Var1[count$Freq==1],]
left=setdiff(TF,data1$TF)
data2=NULL
for(i in 1:length(left)){
  source_temp=data$source[data$TF %in% left[i]]
  source_temp=paste(source_temp,sep="",collapse = ";")
  temp=data.frame(TF=left[i],source=source_temp)
  data2=rbind(data2,temp)
}

result=rbind(data1,data2)  
write.table(result,"TF_union.xls",sep="\t",quote=FALSE,row.names = F)
##---expressed TF----
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/")
TF=read.table("TF/TF_union.xls",sep="\t",header=T,stringsAsFactors = F,quote="")
TCGA=read.table("mRNA/TCGA-STAD-mRNA-tpm.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
TF1=intersect(TCGA$GeneID,TF$TF)
TF=TF[TF$TF %in% TF1,] ##3129
write.table(TF,"TF/TF_union_mRNA.xls",sep="\t",quote=FALSE,row.names = F )
##----TF upset-------
rm(list=ls())
library(UpSetR)         #Upset图（upset 包，适用样本数 2-7）
library(VennDiagram) 
library(ComplexHeatmap)
TF_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/TF"
TF=read.table(file.path(TF_path,"TF_union_mRNA.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/TF")
TFDB=read.table("AnimalATFDB3.0.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
TFDB=unique(TFDB)
colnames(TFDB)=c("TF","source")
hTFtarget=read.table("hTFtarget.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
hTFtarget=unique(hTFtarget)
colnames(hTFtarget)=c("TF","source")
Wei=read.table("Wei.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
Wei$source="Wei"
colnames(Wei)=c("TF","source")
Wei=unique(Wei)
TRRUST=read.table("TRRUST.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
TRRUST$source="TRRUST"
TRRUST=TRRUST[,-2]
colnames(TRRUST)=c("TF","source")
TRRUST=unique(TRRUST)
RcisTarget=read.table("TF_RcisTarget.xls",sep="\t",header=T,stringsAsFactors = F,quote="")
colnames(RcisTarget)=c("TF","source")
RcisTarget=unique(RcisTarget)
hgnc=read.table("hs_hgnc_tfs.txt",sep="\t",header=F,stringsAsFactors = F,quote="")
hgnc$source="hgnc"
colnames(hgnc)=c("TF","source")
hgnc=unique(hgnc)
lt=list(hTFtarget=hTFtarget$TF,
        Wei=Wei$TF,
        TRRUST=TRRUST$TF,
        RcisTarget=RcisTarget$TF,
        TFDB=TFDB$TF,
        HGNC=hgnc$TF)
# aa=Reduce(intersect,list(hTFtarget=hTFtarget$TF,
#                       TRRUST=TRRUST$TF,
#                       RcisTarget=RcisTarget$TF,
#                       TFDB=TFDB$TF))

p=upset(fromList(lt), order.by = "freq",
        nsets = 6,##几个集合
        # nintersects=6,##几个交集
        mb.ratio = c(0.55, 0.45)## 用于控制上下图形所占比例
)
setwd(TF_path)
pdf("TF_upset.pdf",width=6,height=5)
p
dev.off()

###################
##----三、获得STAD 甲基化数据------
##----TCGA STAD hg38版本的甲基化数据--------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/meth")
library(TCGAbiolinks)
##--hg38版本的DNA甲基化数据---#
query_met.hg38 <- GDCquery(project="TCGA-STAD",
                           data.category="DNA Methylation",
                           data.type = "Methylation Beta Value",
                           platform="Illumina Human Methylation 450")
GDCdownload(query_met.hg38)
data.hg38 <- GDCprepare(query_met.hg38)
library(SummarizedExperiment)
hg38=SummarizedExperiment::assay(data.hg38)
saveRDS(data.hg38,"TCGA-STAD_450_hg38.rds")
##############################
##---TCGA STAD hg38版本的甲基化 sample type------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/meth")
TCGA=readRDS("TCGA-STAD_450_hg38.rds")#
TCGA=SummarizedExperiment::assay(TCGA)
TCGA=as.data.frame(TCGA)
group=sapply(strsplit(colnames(TCGA),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
samplesNT=colnames(TCGA)[which(group==1)]
samplesTP=colnames(TCGA)[which(group==0)]
type=data.frame(sample=c(samplesNT,samplesTP),type=c(rep("Normal",length(samplesNT)),
                                                     rep("Tumor",length(samplesTP))))
write.table(type,"sample_type_TCGA.txt",sep="\t",quote=FALSE,row.names = F)
#############################
##-----获得hg38的启动子区域------------
rm(list=ls())
library(ELMER)
library(sesameData)
library(DT)
distal.probes <- get.feature.probe(
  genome = "hg38", # hg38 (default) or hg19
  met.platform = "450K", #  EPIC or 450K (default)
  TSS.range = list(upstream = 2000, downstream = 500), # define promoter regions定义启动子区域
  promoter = TRUE, # 是否选择启动子
  rm.chr = paste0("chr",c("X","Y","M"))
) # 移除的染色体

# #load("C:/Users/admin/AppData/Local/R/win-library/4.2/ELMER.data/data/hm450.hg19.manifest.rda")
# load("C:/Users/admin/AppData/Local/R/win-library/4.2/ELMER.data/data/hm450.hg38.manifest.rda")
# gene=hm450.hg38.manifest@elementMetadata@listData$gene
# chr=hm450.hg38.manifest@elementMetadata@listData$chrm_A
# snp=hm450.hg38.manifest@elementMetadata@listData$MASK_snp5_common
# cg=hm450.hg38.manifest@ranges@NAMES
# chromStart=hm450.hg38.manifest@ranges@start
# chromEnd=hm450.hg38.manifest@ranges@start+hm450.hg38.manifest@ranges@width-1
# data=data.frame(gene=gene,cg=cg,chr=chr,chromStart=chromStart,chromEnd=chromEnd,snp=snp)
# de_snp=data[!data$snp,]
gene=distal.probes@elementMetadata@listData$gene
chr=distal.probes@elementMetadata@listData$chrm_A
snp=distal.probes@elementMetadata@listData$MASK_snp5_common
cg=distal.probes@ranges@NAMES
start=distal.probes@ranges@start
data=data.frame(gene=gene,cg=cg,chr=chr,start=start,snp=snp)
data=data[!data$snp,]##去掉snp位点
loc=grep(";",data$gene)##去掉一对多的甲基化位点
promoter=data[-loc,]
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/meth")
write.table(promoter,"promoter_cg_hg38.xls",sep="\t",quote=FALSE,row.names = F)
# hg38=read.table("illuminaMethyl450_hg38_GDC",sep="\t",header=F,stringsAsFactors = F,quote="")
# colnames(hg38)[1]="cg"
# data2=merge(hg38,promoter,by="cg")
# data3=merge(hg38,data,by="cg")
##############################
##----combat 处理GEO和TCGA甲基化数据的批次效应和knn补缺失值-------
rm(list=ls())
library(sva)
library(DMwR2)
##---获得启动子区域的数据---#
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/meth")
GSE40279=read.table("GSE40279_average_beta.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
TCGA=readRDS("TCGA-STAD_450_hg38.rds")#
TCGA=SummarizedExperiment::assay(TCGA)
TCGA=as.data.frame(TCGA)
del_row=which(rowSums(is.na(TCGA))/ncol(TCGA) > 0.5)
TCGA=TCGA[-del_row, ]
promoter=read.table("promoter_cg_hg38.xls",sep="\t",header=T,quote="")
inter=Reduce(intersect,list(promoter$cg,rownames(TCGA),GSE40279$ID_REF))
TCGA=TCGA[inter,]
rownames(GSE40279)=GSE40279$ID_REF
GSE40279=GSE40279[,-1]
GSE40279=GSE40279[inter,]
TCGA_sample_type=read.table("sample_type_TCGA.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
GSE40279_sample_type=data.frame(sample=colnames(GSE40279),type="Normal")
TCGA=TCGA[,TCGA_sample_type$sample]
TCGA_sample_type$batch=1
GSE40279_sample_type$batch=2
pheno=rbind(TCGA_sample_type,GSE40279_sample_type)
table(pheno$type)##658 normal 395Tumor
saveRDS(GSE40279,"GSE40279_promoter.rds")
##---合并数据 knn补缺失值----#
edata=cbind(TCGA,GSE40279)
pheno$type=factor(pheno$type,levels=c("Tumor","Normal"))
pheno=pheno[order(pheno$batch,pheno$type),]
edata=edata[,pheno$sample]
data1=knnImputation(edata)
#saveRDS(data1,"TCGA_GSE40279_promoter_knn.rds")
anyNA(data1)
##---开始处理批次效应---#
edata=data1
pheno1=pheno
rownames(pheno1)=pheno1$sample
pheno1$sample=1:dim(pheno1)[1]
# 提取batch的信息
batch = pheno$batch
# 建立对照组模型，这不也可以省略
modcombat = model.matrix(~as.factor(type), data=pheno1)
# 去除批次效应
combat_edata1 = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE,  prior.plots=FALSE)
combat_edata1=as.data.frame(combat_edata1)
saveRDS(combat_edata1,"TCGA_GSE40279_promoter_knn_combat.rds")
write.table(pheno,"sample_type.txt",sep="\t",quote=FALSE,row.names = F)
#####################
##-----按照8比2的比例划分甲基化训练集和测试集样本----------
rm(list=ls())
library(caret)
set.seed(8888)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/meth")
data=readRDS("TCGA_GSE40279_promoter_knn_combat.rds")
sample_type=read.table("sample_type.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
Normal=data[,sample_type$sample[sample_type$type=="Normal"]]
Tumor=data[,sample_type$sample[sample_type$type=="Tumor"]]
data1=as.data.frame(t(data))
data1=cbind(sample=rownames(data1),data1)
data1=merge(sample_type,data1,by="sample")
rownames(data1)=data1$sample
data1=data1[,-c(1,3)]
intrain=createDataPartition(y=data1$type,p = 0.8, list = FALSE)
data1=data1[,-1]
train=data1[intrain,]
test=data1[-intrain,]
train=as.data.frame(t(train))
test=as.data.frame(t(test))
# write.table(train,"train.txt",sep="\t",quote=FALSE,row.names = T)
# write.table(test,"test.txt",sep="\t",quote=FALSE,row.names = T)
train_sample=sample_type[sample_type$sample  %in% colnames(train),]
test_sample=sample_type[sample_type$sample %in% colnames(test),]
write.table(train_sample,"train_sample.txt",sep="\t",quote=FALSE,row.names = F)
write.table(test_sample,"test_sample.txt",sep="\t",quote=FALSE,row.names = F)
saveRDS(test,"test_data.rds")
saveRDS(train,"train_data.rds")
#######################
##-----四、获得差异的mRNA和甲基化位点--------
##------获得差异的mRNA------
rm(list=ls())
###--edgeR function------
edgeR<-function(rowcount,condition,comparename,fc,p,q,i,design){
  library(edgeR)
  library(argparser)
  library(base)
  
  compares <- strsplit(comparename,'vs')[[1]]
  condition <-cbind(condition[,1],condition[,i])
  condition <-as.data.frame(condition)
  colnames(condition)<-c("sample","groups1")
  
  if ( colnames(condition)[ncol(condition)]=='type' ) { ngroup <- length(colnames(condition))-2}
  if ( colnames(condition)[ncol(condition)]!='type' ) { ngroup <- length(colnames(condition))-1}
  if ( ngroup==1 ){
    groupdata1 <- condition[which(condition$groups1 %in% compares),]
    colnames(groupdata1)[2] <- 'groups'
    groupdata <- groupdata1 }
  
  if ( ngroup==2 ) {
    groupdata1 <- condition[which(condition$groups1 %in% compares),-3]
    groupdata2 <- condition[which(condition$groups2 %in% compares),-2]
    colnames(groupdata1)[2] <- 'groups'
    colnames(groupdata2)[2] <- 'groups'
    groupdata <- rbind(groupdata1,groupdata2)}
  
  
  rownames(groupdata) <- groupdata[,1]
  group <- factor(groupdata$groups)
  
  if ( design=='normal') {
    designs <- model.matrix(~group)}
  if ( design=='pair' ) {
    type <- factor(groupdata$type)
    designs <- model.matrix(~type+group)}
  
  counts <- subset(rawcount, select=rownames(groupdata))
  counts <- na.omit(counts)
  counts <- round(counts)
  counts <- counts[rowSums(counts)>1,]
  
  y <- DGEList(counts=counts, group=group)
  
  keep <- rowSums(cpm(y) > 1 ) >= 2
  y <- y[keep, ,keep.lib.sizes = FALSE]
  
  # keep <- filterByExpr(y)##filter low reads
  # table(keep)
  # y <- y[keep, , keep.lib.sizes=FALSE]
  
  
  y <- calcNormFactors(y)##TMM normalization is applied to this dataset to account for compositional difference between the libraries.
  #plotMDS(y)
  rownames(designs) <- colnames(y)
  
  if ( design=='normal' ) {
    
    if ( table(groupdata$group)[[compares[1]]] + table(groupdata$group)[[compares[2]]] > 2 ) {
      y <- estimateCommonDisp(y)
      y <- estimateTagwiseDisp(y)
      exact <- exactTest(y, pair=c(compares[2],compares[1]))
      result <- topTags(exact,n=NULL,sort.by='none')}
    
    if ( table(groupdata$group)[[compares[1]]] + table(groupdata$group)[[compares[2]]] == 2 ) {
      exact <- exactTest(y, pair=c(compares[2],compares[1]),dispersion=0.04)
      result <- topTags(exact,n=NULL,sort.by='none')
      y <- equalizeLibSizes(y)}}
  if ( design=='pair') {
    y <- equalizeLibSizes(y)
    y <- estimateGLMCommonDisp(y,designs)
    y <- estimateGLMTrendedDisp(y,designs)
    y <- estimateGLMTagwiseDisp(y,designs)
    
    fit <- glmFit(y, designs)
    lrt <- glmLRT(fit)
    result <- topTags(lrt,n=NULL,sort.by='none') }
  
  
  ID <- rownames(result)
  Count <- y$pseudo.counts
  Count <- Count[ID,]
  
  result <- cbind(ID,as.data.frame(Count),as.data.frame(result))
  
  
  if ( design=='normal' ) {
    result <- subset(result, select=-logCPM)
  }
  if ( design=='pair' ) {
    result <- subset(result,select=-c(logCPM,LR))}
  
  names(result)[names(result)=="logFC"] <- 'log2FoldChange'
  names(result)[names(result)=="PValue"] <- 'pvalue'
  names(result)[names(result)=="FDR"] <- 'padj'
  result$padj[is.na(result$padj)]  <- 1
  result <- result[order(result$pvalue),]
  
  if ( p==0 ) {
    ALL <- subset(result,padj < q & abs(log2FoldChange) > log(fc,2))
    UP <- subset(ALL,log2FoldChange > log(fc,2))
    DOWN <- subset(ALL,log2FoldChange < -log(fc,2))}
  if ( q==0 ) {
    ALL <- subset(result,pvalue <= p & abs(log2FoldChange) >= log(fc,2))
    UP <- subset(ALL,log2FoldChange > log(fc,2))
    DOWN <- subset(ALL,log2FoldChange < -log(fc,2))}
  dir.create("edgeR")
  write.table(result,file=paste("edgeR/",comparename,'_ALL.xls',sep=''),sep='\t',quote=F,row.names=F)
  write.table(ALL,file=paste("edgeR/",comparename,'_DEGs.xls',sep=''),sep='\t',quote=F,row.names=F)
  write.table(UP,file=paste("edgeR/",comparename,'_Up_DEGs.xls',sep=''),sep='\t',quote=F,row.names=F)
  write.table(DOWN,file=paste("edgeR/",comparename,'_Down_DEGs.xls',sep=''),sep='\t',quote=F,row.names=F)
}
#########################
##----差异的mRNA------
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA"
dir.create("D:/project/07.STAD_subtype/analysis0507_v7/01.mRNA_diff")
setwd("D:/project/07.STAD_subtype/analysis0507_v7/01.mRNA_diff")
###class为标签，样本顺序与rawcount一致
condition=read.table(file.path(data_path,"sample_type.txt"), header=T,stringsAsFactors = F,sep="\t")
rawcount=read.table(file.path(data_path,"TCGA-STAD-mRNA-counts.txt"), header=T, row.names=1,sep="\t",check.names = F,quote="")
sample=condition$sample
rawcount=rawcount[,sample]
aa=apply(rawcount,1,sum)
loc=which(aa==0)
rawcount=rawcount[-loc,]
bb=apply(rawcount,1,function(x){loc=which(x==0);len=length(loc);return(len)})
rawcount=rawcount[-which(bb>dim(rawcount)[2]*0.5),]
fc=2
p=0
q=0.05
design="normal"
###class的第i列
i=2
comparename<-'TumorvsNormal'
edgeR(rawcount,condition,comparename,fc,p,q,i,design)
#################
####---绘制火山图------
rm(list=ls())
library(ggplot2)
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA"
setwd("D:/project/07.STAD_subtype/analysis0507_v7/01.mRNA_diff/edgeR")
data=read.delim("TumorvsNormal_ALL.xls",sep="\t",header=T,stringsAsFactors = F,check.names = F)
rownames(data)=data$ID
data=data[,-1]
data$sign="Stable"
up_row=which(data$log2FoldChange>1&data$padj<0.05)
up_mRNA=rownames(data)[up_row]
data[up_row,dim(data)[2]]="Up"#
down_row=which(data$log2FoldChange<c(-1)&data$padj<0.05)
down_mRNA=rownames(data)[down_row]#
data[down_row,dim(data)[2]]="Down" #
limit=max(abs(data$log2FoldChange))+1
my_theme = function(){
  theme(panel.grid = element_blank(),###
        panel.grid.major = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12,hjust=1),
        axis.title = element_text(size=12),##
        plot.title = element_text(hjust=0.5,size=12),#
        plot.margin=unit(c(1,1,1,1), 'lines'))##
}
pdf("huoshantu.pdf",height=4.5,width=5)
#par(mar=c(3,3,1,1),xpd=TRUE)
ggplot(data,aes(x=log2FoldChange,y=-log10(padj),color=sign))+geom_point(size=0.6)+
  scale_color_manual(values=c("#0072B5FF","SeaGreen","#BC3C29FF"))+
  labs(title = "mRNA",colour = "type")+ theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+my_theme()+
  geom_hline(aes(yintercept=-log10(0.05)),linetype="dashed") + 
  geom_vline(aes(xintercept=log(2,2)),linetype="dashed")+ 
  geom_vline(aes(xintercept=log(1/2,2)), linetype="dashed")
dev.off()
###################
##---绘制差异的mRNA的热图----
rm(list=ls())
library(magrittr)
library(ComplexHeatmap)
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA"
setwd("D:/project/07.STAD_subtype/analysis0507_v7/01.mRNA_diff/edgeR")
data=read.delim("TumorvsNormal_DEGs.xls",sep="\t",header=T,stringsAsFactors=F,check.names=F)
row.names(data)=data[,1]
data=data[,-1]
condition=read.table(file.path(data_path,"sample_type.txt"), header=T,stringsAsFactors = F,sep="\t")
table(condition$type)
Up=data%>%dplyr::filter(log2FoldChange > 1)
Down=data%>%dplyr::filter(log2FoldChange < -1)
expr=rbind(Up,Down)
expr=expr[,condition$sample]
###---标准化----###
scaled_expr=t(scale(t(expr)))
ha1=HeatmapAnnotation(df=data.frame(
  Samples=condition$type),
  col=list(Samples=c("Tumor"= "tomato", "Normal"="steelblue")))
ha2=rowAnnotation(df=data.frame(
  mRNA=c(rep("Up",nrow(Up)),rep("Down",nrow(Down)))),
  col=list(mRNA=c("Up"="tomato","Down"="steelblue")))
library(circlize)
pdf(file="complexheatmap_nocluster.pdf",width=4.5,height=3)
Heatmap(scaled_expr,name="expression",top_annotation=ha1,left_annotation =ha2 ,
        cluster_rows=FALSE,cluster_columns=FALSE,show_row_names=FALSE,
        show_column_names=FALSE,
        column_names_gp=gpar(fontsize=12))
dev.off()
pdf(file="complexheatmap_cluster.pdf",width=4.5,height=3)
Heatmap(scaled_expr,name="expression",top_annotation=ha1,left_annotation =ha2 ,
        cluster_rows=FALSE,cluster_columns=TRUE,show_row_names=FALSE,
        show_column_names=FALSE,
        column_names_gp=gpar(fontsize=12))
dev.off()

# bk <- c(seq(-2,4,by=0.001)) #定义热图中值的区间
# pdf(file="complexheatmap_cluster.pdf",width=5,height=3)
# Heatmap(scaled_expr,name="Color_key",top_annotation=ha1,left_annotation =ha2 ,
#         cluster_rows=FALSE,
#         #cluster_columns=FALSE,
#         color=colorRampPalette(rev(c("#D64C49","#FCA86B","#FCFCE3","#8EBCD9","#4979B6")))(length(bk)),
#         show_row_names=FALSE,
#         show_column_names=FALSE,
#         column_names_gp=gpar(fontsize=7))
# dev.off()

# 
# library(pheatmap)
# annotation_col = data.frame( Samples = condition$type)
# rownames(annotation_col) = condition$sample
# annotation_row = data.frame( Genes=c(rep("Up",nrow(Up)),rep("Down",nrow(Down))))
# rownames(annotation_row) = c(rownames(Up),rownames(Down))
# ann_colors = list( Samples = c(Tumor ="#0072B5FF", Normal = "#BC3C29FF"),
#                    Genes=c(Up="#0072B5FF",Down="#BC3C29FF"))
# bk <- c(seq(-2,20,by=0.001)) #定义热图中值的区间
# pdf(file="pheatmap_nocluster.pdf",width=8,height=5)
# pheatmap(scaled_expr,scale = "row",cluster_row = F, cluster_col = F, border=NA,
#          #display_numbers = pmt,
#          show_colnames = F,
#          show_rownames = F,
#          color = colorRampPalette(rev(c("#D64C49","#FCA86B","#FCFCE3","#8EBCD9","#4979B6")))(length(bk)),
#          fontsize_number =8, number_color = "white",
#          annotation_row = annotation_row,
#          annotation_col =annotation_col, 
#          annotation_colors =ann_colors)
# dev.off()
##########################
##=----差异和TF交集------
rm(list=ls())
library(VennDiagram)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/01.mRNA_diff/")
edgeR=read.delim("edgeR/TumorvsNormal_DEGs.xls",sep="\t",header=T,stringsAsFactors=F,check.names=F)
Up=edgeR$ID[edgeR$log2FoldChange>1&edgeR$padj<0.05]
Down=edgeR$ID[edgeR$log2FoldChange<c(-1)&edgeR$padj<0.05]
TF_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/TF"
TF=read.table(file.path(TF_path,"TF_union_mRNA.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
venn_list <- list(Down = Down,TF=TF$TF,Up=Up)
venn.plot=venn.diagram(venn_list, filename = NULL, 
                       fill = c('#0072B5FF', "SeaGreen",'#BC3C29FF'),alpha = 0.50, cat.col = rep('black', 3), 
                       col = 'black', cex = 1, fontfamily = 'serif', 
                       cat.cex = 1, cat.fontfamily = 'serif', scaled = FALSE)

dir.create("DE_TF")
setwd("DE_TF")
pdf(file="venn_DE_TF.pdf",width=5,height=5)
par(mar=c(1,3, 1,3), xpd=TRUE)
grid.draw(venn.plot)
dev.off()
inter=intersect(edgeR$ID,TF$TF)
TF=TF[TF$TF %in% inter,]
write.table(TF,"DE_TF.xls",sep="\t",quote=FALSE,row.names = F)
##########################
##---获得两两差异的甲基化位点 t&chazhi------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/meth")
data=readRDS("train_data.rds")
sample_type=read.table("train_sample.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
Normal=data[,sample_type$sample[sample_type$type=="Normal"]]
Tumor=data[,sample_type$sample[sample_type$type=="Tumor"]]
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
TumorvsNormal=get_diff_point(Normal,Tumor)
sig_T_N=rownames(TumorvsNormal)[which(abs(TumorvsNormal$chazhi)>=0.2 & TumorvsNormal$padj<0.05)]#4716
sig_T_N2=rownames(TumorvsNormal)[which(abs(TumorvsNormal$chazhi)>=0.1 & TumorvsNormal$padj<0.05)]#
sig=TumorvsNormal[sig_T_N,]###4796
fc_chaizhi=TumorvsNormal[,c((dim(TumorvsNormal)[2]-2) : dim(TumorvsNormal)[2])]
fc_chaizhi=cbind(ID=rownames(fc_chaizhi),fc_chaizhi)
diff_meth=list(TumorvsNormal=TumorvsNormal,sig_T_N=sig_T_N,sig=sig,fc_chaizhi=fc_chaizhi)
dir.create("D:/project/07.STAD_subtype/analysis0507_v7/02.meth_diff")
setwd("D:/project/07.STAD_subtype/analysis0507_v7/02.meth_diff")
saveRDS(TumorvsNormal,"TumorvsNormal.rds")
saveRDS(diff_meth,"diff_meth_result.rds")
write.table(sig,"sig_meth.xls",sep="\t",quote=FALSE,row.names=T)
write.table(fc_chaizhi,"fc_chazhi.xls",sep="\t",quote=FALSE,row.names=F)
##---看差异的甲基化位点都是哪些基因---#
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/02.meth_diff")
sig=read.table("sig_meth.xls",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
platform_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/meth"
id=read.table(file.path(platform_path,"promoter_cg_hg38.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
sig_gene=id[id$cg %in% rownames(sig),]
gene=unique(sig_gene$gene)##2311
aa=as.data.frame(table(sig_gene$gene))
write.table(gene,"sig_cg_gene.xls",sep="\t",quote=FALSE,row.names=F)
write.table(sig_gene,"sig_cg_gene_point.xls",sep="\t",quote=FALSE,row.names=F)

##----差异甲基化位点热图--------
rm(list=ls())
library(magrittr)
library(ComplexHeatmap)
data_path="D:/project/07.STAD_subtype/analysis0507_v7//00.data/meth"
setwd("D:/project/07.STAD_subtype/analysis0507_v7/02.meth_diff")
data=readRDS("diff_meth_result.rds")
TumorvsNormal=data$TumorvsNormal
sig=data$sig
condition=read.table(file.path(data_path,"train_sample.txt"), header=T,stringsAsFactors = F,sep="\t")
condition=condition[condition$sample %in% colnames(TumorvsNormal),]
Up=sig%>%dplyr::filter(chazhi >=0.2)
Down=sig%>%dplyr::filter(chazhi <= -0.2)
expr=rbind(Up,Down)
expr=expr[,condition$sample]
###---标准化----###
scaled_expr=t(scale(t(expr),
                    center = T, scale=T)) 
ha1=HeatmapAnnotation(df=data.frame(
  Samples=condition$type),
  col=list(Samples=c("Tumor"= "tomato", "Normal"="steelblue")))
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

# bk <- c(seq(-2,4,by=0.001)) #定义热图中值的区间
# pdf(file="complexheatmap_cluster.pdf",width=5,height=3)
# Heatmap(scaled_expr,name="Color_key",top_annotation=ha1,left_annotation =ha2 ,
#         cluster_rows=FALSE,
#         #cluster_columns=FALSE,
#         color=colorRampPalette(rev(c("#D64C49","#FCA86B","#FCFCE3","#8EBCD9","#4979B6")))(length(bk)),
#         show_row_names=FALSE,
#         show_column_names=FALSE,
#         column_names_gp=gpar(fontsize=7))
# dev.off()

# 
library(pheatmap)
annotation_col = data.frame( Samples = condition$type)
rownames(annotation_col) = condition$sample
annotation_row = data.frame( CpG=c(rep("Up",nrow(Up)),rep("Down",nrow(Down))))
rownames(annotation_row) = c(rownames(Up),rownames(Down))
ann_colors = list( Samples = c(Tumor ="#BC3C29FF", Normal = "#0072B5FF"),
                   CpG=c(Up="#BC3C29FF",Down="#0072B5FF"))
scaled_expr[scaled_expr>5]=5
scaled_expr[scaled_expr< -5]= -5
bk <- c(seq(-5,5,by=0.1)) #定义热图中值的区间
pdf(file="pheatmap_cluster.pdf",width=5,height=5)
pheatmap(scaled_expr,
         #scale = "row",
         cluster_row = F, cluster_col = T, border=NA,
         #display_numbers = pmt,
         show_colnames = F,
         show_rownames = F,
         color = colorRampPalette(rev(c("#BC3C29FF","#FCFCE3","#0072B5FF")))(length(bk)),
         fontsize_number =8, number_color = "white",
         annotation_row = annotation_row,
         annotation_col =annotation_col, 
         annotation_colors =ann_colors)
dev.off()
################