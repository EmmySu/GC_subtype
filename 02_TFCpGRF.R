##---六、WGCNA获得关键模块和重要的转录因子-------
##----所有RNA-seq WGCNA------
rm(list=ls())
library(WGCNA)
set.seed(666)
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA"
expression=read.table(file.path(data_path,"TCGA-STAD-mRNA-counts.txt"),sep="\t",header=T,row.names=1, check.names=F,quote="")
sample_type=read.table(file.path(data_path,"sample_type.txt"),sep="\t",header=T, check.names=F,quote="")
sample_type=sample_type[sample_type$type=="Tumor",]
expression=expression[,sample_type$sample]
aa=apply(expression,1,sum)
loc=which(aa==0)
expression=expression[-loc,]
bb=apply(expression,1,function(x){loc=which(x==0);len=length(loc);return(len)})
expression=expression[-which(bb>dim(expression)[2]*0.5),]
# expression=expression[rowMeans(expression)>0.1,]
# expression=expression[apply(expression,1,sd)>0.1,] #删除波动小的基因
#Step1 对表达数据进行预处理，去掉不合格的点，根据层次聚类，得到离群数据，去掉离群样本,获得最终的表达数据为datExpr
datExpr0 = as.data.frame(t(expression))
expression=cbind(ID=rownames(expression),expression)
###使用goodSamplesGenes筛选gene##
#对datExpr0的数据进行goodSamplesGenes测试，看符合标准的基因的返回值
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

#如果有不合格的基因或数据，则显示出来，并剔除
if (!gsg$allOK)
{ # 显示出要被删除的基因和数据
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # 去掉不合格的基因和数据
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
##
dir.create("D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count")
setwd("D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count")
##样品聚类,查找游离的样品
sampleTree = hclust(dist(datExpr0), method = "average")
###
pdf(file = "sampleClustering.pdf", width = 60, height = 5);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()
###删除剪切线以下的样品
clust=cutreeStatic(sampleTree,cutHeight = 4e+06, minSize = 10)
table(clust)
keepSamples=(clust==1)
##样本都保留，没有离群样本
datExpr0 = datExpr0[keepSamples,]
write.table(datExpr0, file = "datExpr0.txt",sep="\t",row.names=T,quote=F)
##Step2 选择软阈值##
#设定画图显示时的阈值，1~20
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# 调用网络拓扑分析功能
sft <- pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)##在sft这个对象中保存了每个power值计算出来的网络的特征
#其中sft$powerEstimate就是最佳的power值，sft$fitIndices保存了每个power对应的网络的特征。我们可以对不同power取值下，网络的特征进行可视化，代码如下
#对结果作图
cex1 <- 0.9##这个值是可改的，一般要求k与p(k)的相关性达到0.85时的power作为β值,这里取0.9
pdf("soft_Threshold.pdf",height=10,width=20)
par(mfrow = c(1,2))
# 计算软阈值，并画图展示，设定abline寻找最适合的阈值
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 根据R^2值为0.9时 筛选阈值
abline(h = 0.9,col="red")##
# 同时观察不同阈值的连接度,sft$fitIndices保存了每个power构建的相关性网络中的连接度的统计值，k就是连接度值，这里对连接度的均值进行可视化
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()
##邻接矩阵转换
sft #查看最佳power值
softPower =sft$powerEstimate #最佳power值
adjacency = adjacency(datExpr0, power = softPower)
softPower
###TOM矩阵
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
###基因聚类
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file="gene_clustering.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()
###动态剪切模块识别
minModuleSize=60      #模块基因数目
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

#绘制基因模块的图形
pdf(file="Dynamic_Tree.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
###查找相似模块聚类
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file="Clustering_module.pdf",width=7,height=7)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25    #模块相似性的阈值
abline(h=MEDissThres, col = "red")
dev.off()
###相似模块合并
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
#合并模块的图形
pdf(file="merged_dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
###基因所在的模块
probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "module_all.txt",sep="\t",row.names=F,quote=F)
###输出每个模块的基因
mm=unique(moduleColors)

for(i in 1:length(mm))
{
  modules <-mm[i]##可改
  # 选择模块中基因
  sybmbol = colnames(datExpr0)
  inModule = is.finite(match(moduleColors, modules)) #isFinite() 函数用于检查其参数是否是无穷大
  modProbes = sybmbol[inModule]
  # modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  
  # Export the network into edge and node list files Cytoscape can read
  # threshold 默认为0.5, 可以根据自己的需要调整，也可以都导出后在 cytoscape中再调整
  dir.create("modules")
  dir.create(paste0("modules/",modules,sep=""))
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("modules/",modules,"/edges.txt",sep=""),
                                 nodeFile = paste("modules/",modules,"/nodes.txt",sep=""),
                                 
                                 threshold = 0.2,
                                 nodeNames = modProbes,
                                 nodeAttr = moduleColors[inModule])
  
}

mm=unique(moduleColors)

for(i in 1:length(mm))
{
  modules <-mm[i]##可改
  # 选择模块中基因
  sybmbol = colnames(datExpr0)
  inModule = is.finite(match(moduleColors, modules)) #isFinite() 函数用于检查其参数是否是无穷大
  modProbes = sybmbol[inModule]
  # modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  # threshold 默认为0.5, 可以根据自己的需要调整，也可以都导出后在 cytoscape中再调整
  dir.create("modules")
  dir.create(paste0("modules/",modules,sep=""))
  write.table(modProbes,paste0("modules/",modules,"/module_nodes.txt",sep=""),sep="\t",quote=FALSE,row.names = F)
  
}
###基因所在的模块
probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "module_all.txt",sep="\t",row.names=F,quote=F)
###输出每个模块的基因
for (mod in 1:nrow(table(moduleColors))){  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("module_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}
save(expression,datExpr0,moduleColors,modules,sft,TOM ,file = "WGCNA_result.RData")##保存数据
##################
#################
##----挑选重要的模块-----------
##---查看不同模块有多少转录因子 甲基化位点对应的基因 有多少个是显著的 ------
rm(list=ls())
library(cowplot)
library(VennDiagram)
TF_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/TF"
TF_sig_path="D:/project/07.STAD_subtype/analysis0507_v7/01.mRNA_diff/DE_TF"
cg_path="D:/project/07.STAD_subtype/analysis0507_v7/02.meth_diff"
TF=read.table(file.path(TF_path,"TF_union_mRNA.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
TF_sig=read.table(file.path(TF_sig_path,"DE_TF.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
cg_gene=read.table(file.path(cg_path,"sig_cg_gene.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
colnames(cg_gene)[1]="gene"
cg_gene=cg_gene[!is.na(cg_gene$gene),]
cg_gene=as.data.frame(cg_gene)
colnames(cg_gene)[1]="gene"
setwd("D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count/")
node_all=read.table("module_all.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
module=unique(node_all$moduleColor)
result=NULL
result2=NULL
dir.create("TF_cg_inter")
plot_list=list()
plot_list1=list()
for(i in 1:length(module)){
  node=node_all[node_all$moduleColor==module[i],]
  TF_inter=intersect(node$probes,TF$TF)
  TF_inter_sig=intersect(node$probes,TF_sig$TF)
  cg_inter=intersect(node$probes,cg_gene$gene)
  if(length(TF_inter_sig)<0){
    TF_inter_sig=NULL
  }
  if(length(cg_inter)<0){
    cg_inter=NULL
  }
  venn_list <- list(module=node$probes,TF = TF$TF,TF_sig=TF_sig$TF,cg_gene=cg_gene$gene)
  names(venn_list)[1]=module[i]
  venn.plot=venn.diagram(venn_list, filename = NULL, 
                         fill = c('#0072B5FF', '#BC3C29FF',"#b2df8a","#F28E2B"),alpha = 0.50, cat.col = rep('black', 4), 
                         col = 'black', cex = 1, fontfamily = 'serif', 
                         cat.cex = 1, cat.fontfamily = 'serif', scaled = FALSE)
  plot_list[[i]]=venn.plot
  venn_list1 <- list(module=node$probes,TF_sig=TF_sig$TF,cg_gene=cg_gene$gene)
  names(venn_list1)[1]=module[i]
  venn.plot1=venn.diagram(venn_list1, filename = NULL, 
                          fill = c('#0072B5FF', '#BC3C29FF',"#b2df8a"),alpha = 0.50, cat.col = rep('black', 3), 
                          col = 'black', cex = 1, fontfamily = 'serif', 
                          cat.cex = 1, cat.fontfamily = 'serif', scaled = FALSE)
  plot_list1[[i]]=venn.plot1
  temp=data.frame(module=module[i],TF=paste(TF_inter,sep="",collapse = ";"),TF_inter_sig=paste(TF_inter_sig,sep="",collapse = ";"),
                  cg_gene=paste(cg_inter,sep="",collapse = ";"))
  temp2=data.frame(module=module[i],TF_inter=length(TF_inter),TF_inter_sig=length(TF_inter_sig),cg_inter=length(cg_inter),node=dim(node)[1])
  result=rbind(result,temp)
  result2=rbind(result2,temp2)
}
pdf("TF_cg_inter/TF_cg_inter.pdf",height=9.5,width=6)
plot_grid(plotlist=plot_list,ncol=3)
dev.off()
pdf("TF_cg_inter/TF_cg_inter1.pdf",height=9.5,width=6)
plot_grid(plotlist=plot_list1,ncol=3)
dev.off()
write.table(result,"TF_cg_inter/module_TF_cg.xls",sep="\t",quote=FALSE,row.names = F)
write.table(result2,"TF_cg_inter/module_TF_cg_number.xls",sep="\t",quote=FALSE,row.names = F)
####################
##----查看模块都对应哪些甲基化位点和转录因子-----
##----并且做模块的甲基化位点和转录因子的相关性-------
rm(list=ls())
library(psych)#psych包用于计算相关性、p值等信息
library(pheatmap)
library(reshape2)
##----function---#
##--获得模块的甲基化位点和转录因子---#
get_cg_TF <- function(module,tf_cg,cg_gene){
  temp= tf_cg[which(tf_cg$module==module),]
  TF=unlist(strsplit(temp$TF_inter_sig[temp$module==module],";"))
  gene=unlist(strsplit(temp$cg_gene[temp$module==module],";"))
  cg=cg_gene$cg[cg_gene$gene %in% gene]
  TF_cg=list(TF=TF,cg_gene=gene,cg=cg)
  return(TF_cg)
}
##--获得模块的甲基化位点和转录因子的相关性---#
cor_TF_cg <- function(meth,GeneExp,TF_cg){
  cg=TF_cg$cg
  TF=TF_cg$TF
  if(length(cg)>0 & length(TF)>0){
    cg=meth[,cg]
    TF=GeneExp[,TF]
    res=corr.test(cg, TF, method = "pearson",adjust="none")
    #提取相关性、p值
    cmt=res$r
    pmt=res$p
    result=melt(cmt,value.name="cor")
    result$pvalue=as.vector(pmt)
    colnames(result)[1]="cg"
    colnames(result)[2]="TF"
    sig=result[which(result$pvalue<0.05 & result$cor<0 ),]
    cg_freq=as.data.frame(table(sig$cg))
    TFfreq=as.data.frame(table(sig$TF))
    length(which(TFfreq$Freq>median(TFfreq$Freq)))
    cor_sig_p=dim(sig)[1]/dim(result)[1]
    cg_sig=cg_freq$Var1[which(cg_freq$Freq>length(which(TFfreq$Freq>median(TFfreq$Freq))))] 
    result=list(result=result,sig=sig,cor_sig_p=cor_sig_p,cg_freq=cg_freq,TFfreq=TFfreq,cg_sig=cg_sig)
    
  }else{
    result=list()
  }
  return(result)
}
##---获得模块的信息---#
path="D:/project/07.STAD_subtype/analysis0507_v7/02.meth_diff"
cg_gene=read.table(file.path(path,"sig_cg_gene_point.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
cg_tf="D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count/TF_cg_inter"
tf_cg_module=read.table(file.path(cg_tf,"module_TF_cg.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
##---读取基因表达和甲基化表达---#
meth_path="D:/project/07.STAD_subtype/analysis0507_v7/02.meth_diff"
mRNA_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA"
meth=readRDS(file.path(meth_path,"diff_meth_result.rds"))
sig=meth$sig
mRNA=read.table(file.path(mRNA_path,"TCGA-STAD-mRNA-counts.txt"),sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F,row.names = 1)
colnames(sig)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(sig))
colnames(mRNA)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4",colnames(mRNA))
sample=intersect(colnames(sig),colnames(mRNA)) ##271
Meth=sig[,sample]
GeneExp=mRNA[,sample]
meth1=t(Meth)
GeneExp1=t(GeneExp)
module=unique(tf_cg_module$module)
all_result=list()
cor_result=NULL
for(i in 1:length(module)){
  TF_cg=get_cg_TF(module[i],tf_cg_module,cg_gene)
  result=cor_TF_cg(meth1,GeneExp1,TF_cg)
  if(length(result)>0){
    all_result[[i]]=result
    names(all_result)[i]=module[i]
    temp=data.frame(module=module[i],cor_p=result$cor_sig_p)
    cor_result=rbind(cor_result,temp)
  }
  print(i)
}
setwd("D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count")
dir.create("find_important")
setwd("find_important")
saveRDS(all_result,"module_cg_gene_cor.rds")
write.table(cor_result,"cor_p.xls",sep="\t",quote=FALSE,row.names = F)
###################
##-----获得每个模块的差异甲基化和差异转录因子的显著性 超几何检验---------
rm(list=ls())
##---筛选显著的模块---#
setwd("D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count/")
DE_path="D:/project/07.STAD_subtype/analysis0507_v7/01.mRNA_diff/edgeR"
TF_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/TF"
cg_path="D:/project/07.STAD_subtype/analysis0507_v7/02.meth_diff"
TF=read.table(file.path(TF_path,"TF_union.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
DE=read.table(file.path(DE_path,"TumorvsNormal_DEGs.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
TF_sig=TF[TF$TF %in% DE$ID,]
cg_gene=read.table(file.path(cg_path,"sig_cg_gene.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
colnames(cg_gene)[1]="gene"
cg_gene=cg_gene[!is.na(cg_gene$gene),]
cg_gene=as.data.frame(cg_gene)
colnames(cg_gene)[1]="gene"
setwd("D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count/")
node_all=read.table("module_all.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
number=read.table("TF_cg_inter/module_TF_cg_number.xls",sep="\t",header=T,stringsAsFactors = F,quote="")
##---TF是不是显著---#
M=length(intersect(TF_sig$TF,node_all$probes))##601
N=length(node_all$probes)#18264
number$TF_pvalue=0
for(i in 1:dim(number)[1]){
  n=number$node[i]
  m=number$TF_inter_sig[i]
  pvalue=phyper(m,M,(N-M),n,lower.tail = F)
  number$TF_pvalue[i]=round(pvalue,3)
}

##---cg是不是显著---#
M=length(intersect(cg_gene$gene,node_all$probes))##1888
N=length(node_all$probes)
number$cg_pvalue=0
for(i in 1:dim(number)[1]){
  n=number$node[i]
  m=number$cg_inter[i]
  pvalue=phyper(m,M,(N-M),n,lower.tail = F)
  number$cg_pvalue[i]=round(pvalue,3)
}
number$TF_FDR=p.adjust(number$TF_pvalue,method="fdr",n=length(number$TF_pvalue))
number$cg_FDR=p.adjust(number$cg_pvalue,method="fdr",n=length(number$cg_pvalue))

loc=which(number$TF_FDR<0.05 & number$cg_FDR<0.05)
module=number[loc,]
number$TF_module=length(intersect(TF_sig$TF,node_all$probes))
number$cg_module=length(intersect(cg_gene$gene,node_all$probes))
number$WGCNA_number=length(node_all$probes)
write.table(number,"TF_cg_inter/module_p.xls",sep="\t",quote=FALSE,row.names = F)
######################
##---提出一个方法考虑TF,cg_gene在模块的出现的概率以及两者的相关性的概率-----
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count/find_important")
p_module=read.table("D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count/TF_cg_inter/module_p.xls",
                    sep="\t",header=T,stringsAsFactors = F,quote="")
cor_result=read.table("cor_p.xls",sep="\t",header=T,stringsAsFactors = F,quote="")
p_module=p_module[,c("module","TF_pvalue","cg_pvalue","TF_FDR","cg_FDR" )]
cor_result=cor_result[order(cor_result$cor_p,decreasing = T),]
cor_result$rank=1:dim(cor_result)[1]
p_module=merge(p_module,cor_result,by="module")
p_module$G_p_TF=1
p_module$G_p_TF[p_module$TF_pvalue>0.05]=0
p_module$G_p_cg=1
p_module$G_p_cg[p_module$cg_pvalue>0.05]=0
p_module$ComP=p_module$G_p_TF*p_module$G_p_cg*p_module$cor_p
p_module=p_module[order(p_module$TF_pvalue),]
p_module$TF_rank=1:dim(cor_result)[1]
p_module=p_module[order(p_module$cg_pvalue),]
p_module$cg_rank=1:dim(cor_result)[1]
p_module$mean_rank=(p_module$TF_rank+p_module$cg_rank+p_module$rank)/3
p_module=p_module[order(p_module$mean_rank),]
p_module$final_rank=1:dim(cor_result)[1]
p_module$p_final=p_module$TF_pvalue*(1-p_module$cor_p)+p_module$cg_pvalue*(1-p_module$cor_p)
p_module$p_final2=p_module$TF_pvalue*(1-p_module$cor_p)+p_module$cg_pvalue*(1-p_module$cor_p)+(1-p_module$cor_p)
p_module$p_final3=p_module$TF_pvalue*(1-p_module$cor_p)+p_module$cg_pvalue*(1-p_module$cor_p)+(1-p_module$cor_p)
p_module$p_final4=p_module$TF_pvalue*(1-p_module$cor_p)+p_module$cg_pvalue*(1-p_module$cor_p)+(1-p_module$ComP)
write.table(p_module,"module_p_final.xls",sep="\t",quote=FALSE,row.names = F)

################
#--比较ComP的性能，选择最佳阈值------
rm(list=ls())
library(reshape2)
library(tidyr)
library(ggplot2)
library(cowplot)
library(gtable)
library(grid)
#----设置多个y轴-----#
hinvert_title_grob <- function(grob){
  # 交换宽度
  widths <- grob$widths
  grob$widths[1] <- widths[3]
  grob$widths[3] <- widths[1]
  grob$vp[[1]]$layout$widths[1] <- widths[3]
  grob$vp[[1]]$layout$widths[3] <- widths[1]
  
  # 修改对齐 
  grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
  grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
  grob$children[[1]]$x <- unit(1, "npc") - grob$children[[1]]$x
  grob
}

add_another_yaxis <- function(g1, g2, offset = 0) {
  # ============ 1. 主绘图区 ============ #
  # 获取主绘图区域
  pos <- c(subset(g1$layout, name == "panel", select = t:r))
  # 添加图形
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], 
                       pos$t, pos$l, pos$b * ((offset - 2) * 0.00001 + 1), pos$l)
  # ============ 2. 轴标签 ============ #
  index <- which(g2$layout$name == "ylab-l")
  ylab <- g2$grobs[[index]]
  ylab <- hinvert_title_grob(ylab)
  # 添加轴标签
  g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], pos$r)
  g <- gtable_add_grob(g, ylab, pos$t, pos$r + 1, pos$b, pos$r + 1, clip = "off", name = "ylab-r")
  # ============ 3. 轴设置 ============ #
  index <- which(g2$layout$name == "axis-l")
  yaxis <- g2$grobs[[index]]
  # 将 Y 轴线移动到最左边
  yaxis$children[[1]]$x <- unit.c(unit(0, "npc"), unit(0, "npc"))
  # 交换刻度线和刻度标签
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  # 移动刻度线
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, "npc") + unit(3, "pt")
  # 刻度标签位置转换和对齐
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  yaxis$children[[2]] <- ticks
  # 添加轴，unit(3, "mm") 增加轴间距
  g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l] , pos$r)
  g <- gtable_add_grob(g, yaxis, pos$t, pos$r + 1, pos$b, pos$r + 1, clip = "off", name = "axis-r")
  g
}
# 接受可变参数，可添加多个 Y 轴
plot_multi_yaxis <- function(..., right_label_reverse = TRUE) {
  args <- list(...)
  len <- length(args)
  g <- ggplotGrob(args[[1]])
  for (i in len:2) {
    if (right_label_reverse) {
      # 为轴标签添加旋转
      args[[i]] <- args[[i]] + theme(axis.title.y = element_text(angle = 270))
    }
    g2 <- ggplotGrob(args[[i]])
    g <- add_another_yaxis(g, g2, offset = i)
  }
  # 绘制图形
  grid.newpage()
  grid.draw(g)
  g
}
##--找模块重要的甲基化位点---#
setwd("D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count/find_important")
all_result=readRDS("module_cg_gene_cor.rds")
cor_p=read.table("module_p_final.xls",sep="\t",header=T,stringsAsFactors = F,quote="")
library(ggplot2)
aa=seq(0,1,0.01)
threshold_result=NULL
for(i in 1:length(aa)){
  threshold=aa[i]
  ll=length(which(cor_p$p_final4<threshold))
  temp=data.frame(threshold=threshold,number=ll)
  threshold_result=rbind(threshold_result,temp)
}
rank=data.frame(module=c(rep(cor_p$module,3)),rank=c(cor_p$rank,cor_p$TF_rank,cor_p$cg_rank),
                type=c(rep("cor",dim(cor_p)[1]),
                       rep("TF",dim(cor_p)[1]),
                       rep("cg",dim(cor_p)[1])))
ComP=data.frame(module=cor_p$module,ComP=cor_p$p_final4,
                type=c(rep("ComP",dim(cor_p)[1])))
ComP=ComP[order(ComP$ComP,decreasing = F),]
rank$module=factor(rank$module,levels = ComP$module)
p1=ggplot(data = rank,aes(x=module,y=rank,group = type,color=type,shape=type))+
  geom_point(size=3)+
  #geom_line(size=1)+
  xlab("")+#横坐标名称
  ylab("Rank")+#纵坐标名称
  #labs(title = patient[i])+
  #scale_fill_brewer(palette = 'Set1')
  #scale_color_manual(values = cols)+
  #scale_x_discrete(expand = c(0.1,0))+
  theme_bw()+
  theme(legend.title=element_text(size=7) ,legend.position = "bottom", legend.text=element_text(size=7))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA),
        plot.title = element_text(hjust=0.5,size=7),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line.y = element_line(color = "black"),
        axis.text.y = element_text(color = "black"), 
        axis.ticks.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ComP$module=factor(ComP$module,levels = ComP$module)
p2=ggplot(data = ComP,aes(x=module,y=ComP,group = type,color=type,shape=type))+
  geom_point(shape=1,size=3)+
  geom_line(size=1)+
  #geom_hline(aes(yintercept=1), colour="#990000", linetype="dashed")+
  xlab("")+#横坐标名称
  ylab("ComP")+#纵坐标名称
  #labs(title = patient[i])+
  #scale_fill_brewer(palette = 'Set1')
  scale_color_manual(values = "purple")+
  #scale_x_discrete(expand = c(0.1,0))+
  theme_bw()+
  theme(legend.title=element_text(size=7) ,legend.position = "bottom", legend.text=element_text(size=7))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA),
        plot.title = element_text(hjust=0.5,size=7),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line.y = element_line(color = "black"),
        axis.text.y = element_text(color = "black"), 
        axis.ticks.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

p2_legend=get_legend(p2+ theme(legend.position = "bottom"))
p=plot_multi_yaxis(p1, p2)
plot2<- plot_grid(p, p2_legend,nrow=2,rel_heights = c(3,0.05))
pdf("important_compare.pdf",width=6.5,height=5)
plot2
dev.off()
#################
##---找模块重要的TF------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count/find_important")
all_result=readRDS("module_cg_gene_cor.rds")
cg_sig=read.table("module_p_final.xls",sep="\t",header=T,stringsAsFactors = F,quote="")
cg_sig=cg_sig[order(cg_sig$p_final4),]
module=cg_sig$module[cg_sig$p_final4<1][1]
result=NULL
for( i in 1:length(module)){
  temp=all_result[[module[i]]]
  sig=temp$sig
  TF_freq=temp$TFfreq
  cg_freq=temp$cg_freq
  TF_sig=TF_freq$Var1[which(TF_freq$Freq>=round(length(cg_freq$Var1)/2))]
  temp_TF=data.frame(module=module[i],TF=TF_sig)
  result=rbind(result,temp_TF)
}
write.table(result,"important_TF.xls",sep="\t",quote=FALSE,row.names = F)
length(unique(result$TF)) ##108
##################
##--找模块重要的甲基化位点-----
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count/find_important")
all_result=readRDS("module_cg_gene_cor.rds")
cg_sig_final=NULL
TF=read.table("important_TF.xls",sep="\t",header=T,stringsAsFactors = F,quote="")
module=unique(TF$module)
for(i in 1:length(module)){
  temp=all_result[[module[i]]]
  sig=temp$sig
  TF_freq=temp$TFfreq
  cg_freq=temp$cg_freq
  cg_sig=cg_freq$Var1[which(cg_freq$Freq>=round(length(TF_freq$Var1)/2))]
  temp_cg=data.frame(module=module[i],cg=cg_sig,ll=length(cg_sig))
  cg_sig_final=rbind(cg_sig_final,temp_cg)
}
write.table(cg_sig_final,"important_cg.xls",sep="\t",quote=FALSE,row.names = F)
##---六、获得单因素显著的TF和甲基化位点---
##---TF是模块中重要的TF，甲基化位点是模块重要的甲基化位点----
##---得到OS TF------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/00.data/")
##--survival--##
df=read.table("clinical/TCGA-STAD.survival.tsv",sep="\t",header=T,stringsAsFactors = F,quote="")
df=df[grep("-01",df$sample),]
df=df[,c(3,2,4)]
colnames(df)[1]="sample"
clinical=df[,c("sample","OS","OS.time")]
rownames(clinical)=clinical$sample
clinical=clinical[,-1]
##--mRNA--##
TF_path="D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count/find_important"
TF=read.table(file.path(TF_path,"important_TF.xls"),sep="\t",header=T,stringsAsFactors = F,quote = "")[,2]
data=read.table("mRNA/TCGA-STAD-mRNA-tpm_log.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
data=data[TF,]
sample_type=read.table("mRNA/sample_type.txt",sep="\t",quote="",header=T,stringsAsFactors = F)
data=data[,sample_type$sample[sample_type$type=="Tumor"]]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))
data=t(data)
sameSample=intersect(row.names(data),row.names(clinical))##279
data=data[sameSample,]
clinical=clinical[sameSample,]
out=cbind(clinical,data)
out=cbind(id=row.names(out),out)
dir.create("D:/project/07.STAD_subtype/analysis0507_v7/04.OS")
setwd("D:/project/07.STAD_subtype/analysis0507_v7/04.OS")
write.table(out,file="OS_gene.txt",sep="\t",row.names=F,quote=F)
##---重要的TF单因素cox回归分析-------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/04.OS")
data=read.table("OS_gene.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
datass=data[,-1]
colnames(datass)[c(1,2)]=c("status","time")
library(survival)
library(plyr)
BaSurv=Surv(time =datass$time,event=datass$status)
datass$BaSurv=with(datass,BaSurv)
UniCox=function(datass,x){
  type=unique(datass[,which(colnames(datass)==x)])
  FML=as.formula(paste0('BaSurv~',x))
  GCox=coxph(FML,data=datass)
  GSum=summary(GCox)
  HR=round(GSum$coefficients[,2],2)
  PValue=round(GSum$coefficients[,5],3)
  HR.95L=round(GSum$conf.int[,3],2)
  HR.95H=round(GSum$conf.int[,4],2)
  CI=paste0(round(GSum$conf.int[,3:4],2),collapse="-")
  UniCox=data.frame('Characteristics'=x,'Hazard Ratio'=HR,'HR.95L'=HR.95L,'HR.95H'=HR.95H,'CI95'=CI,'p value'=PValue)
  return(UniCox)
}

VarNames=colnames(datass)[-c(1,2,dim(datass)[2])]
result=NULL
for( i in 1:length(VarNames)){
  temp=UniCox(datass,VarNames[i])
  result=rbind(result,temp)
}
result$Characteristics[result$p.value<0.05]
UniVar.sig=result[which(result$p.value<0.05),]#37
dir.create("TF")
write.table(result,"TF/single_cox_result.txt",sep="\t",quote=FALSE,row.names = F)
write.table(UniVar.sig,"TF/single_cox_sig.txt",sep="\t",quote=FALSE,row.names = F)
##---绘制单因素显著的TF的森林图----------
rm(list=ls())
library(survival) 
library(survminer)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/04.OS/TF")
coxFile="single_cox_sig.txt"
#绘制森林图函数##
bioForest <- function(coxFile=null, forestFile=null){
  #读取输入文件
  rt=read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
  rt=rt[order(rt$Hazard.Ratio),]
  gene=rownames(rt)
  hr=rt$"Hazard.Ratio"
  hrLow =rt$"HR.95L"
  hrHigh=rt$"HR.95H"
  Hazard.ratio=rt$CI95
  pVal=ifelse(rt$p.value<0.001, "<0.001", sprintf("%.3f", rt$p.value))
  hrLow=as.numeric(hrLow)
  hrHigh=as.numeric(hrHigh)
  hr=as.numeric(hr)
  #输出图形
  pdf(file=forestFile, width=7.6, height=7.5)
  n=nrow(rt)
  nRow=n+1
  ylim=c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim =c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #绘制右边森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim =c(0,max(as.numeric(hrLow),as.numeric(hrHigh))+1)
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#98c1d9",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor =ifelse(as.numeric(rt$Hazard.Ratio) < 1, "#98c1d9", "#FB8072")
  points(as.numeric(hr), n:1, pch =15, col =boxcolor, cex=1.5)
  axis(1)
  dev.off()
}
bioForest(coxFile=coxFile, forestFile="TF_forest.pdf")
#---得到重要的甲基化位点的OS-----------
rm(list=ls())
clinical_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/clinical/"
df=read.table(file.path(clinical_path,"TCGA-STAD.survival.tsv"),sep="\t",header=T,stringsAsFactors = F,quote="")
df=df[grep("-01",df$sample),]
df=df[,c(3,2,4)]
colnames(df)[1]="sample"
clinical=df[,c("sample","OS","OS.time")]
rownames(clinical)=clinical$sample
clinical=clinical[,-1]
meth_path="D:/project/07.STAD_subtype/analysis0507_v7/02.meth_diff"
meth_path2="D:/project/07.STAD_subtype/analysis0507_v7/00.data/meth"
sample_train=read.table(file.path(meth_path2,"train_sample.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
DE=readRDS(file.path(meth_path,"diff_meth_result.rds"))
sig=DE$sig
meth_important_path="D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count/find_important"
meth_im=read.table(file.path(meth_important_path,"important_cg.xls"),sep="\t",header=T,stringsAsFactors = F,quote = "")[,2]
meth=readRDS(file.path(meth_path2,"TCGA_GSE40279_promoter_knn_combat.rds"))
meth=meth[meth_im,sample_train$sample]
meth=meth[,sample_train$sample[sample_train$type=="Tumor"]]
colnames(meth)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(meth))
meth=t(meth)
sameSample=intersect(row.names(meth),row.names(clinical))##306
meth=meth[sameSample,]
clinical=clinical[sameSample,]
out=cbind(clinical,meth)
out=cbind(id=row.names(out),out)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/04.OS")
write.table(out,file="OS_CpG.txt",sep="\t",row.names=F,quote=F)
##---univar CpG------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/04.OS")
data=read.table("OS_CpG.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
datass=data[,-1]
colnames(datass)[c(1,2)]=c("status","time")
library(survival)
library(plyr)
BaSurv=Surv(time =datass$time,event=datass$status)
datass$BaSurv=with(datass,BaSurv)
UniCox=function(datass,x){
  type=unique(datass[,which(colnames(datass)==x)])
  FML=as.formula(paste0('BaSurv~',x))
  GCox=coxph(FML,data=datass)
  GSum=summary(GCox)
  HR=round(GSum$coefficients[,2],2)
  PValue=round(GSum$coefficients[,5],3)
  HR.95L=round(GSum$conf.int[,3],2)
  HR.95H=round(GSum$conf.int[,4],2)
  CI=paste0(round(GSum$conf.int[,3:4],2),collapse="-")
  UniCox=data.frame('Characteristics'=x,'Hazard Ratio'=HR,'HR.95L'=HR.95L,'HR.95H'=HR.95H,'CI95'=CI,'p value'=PValue)
  return(UniCox)
}

VarNames=colnames(datass)[-c(1,2,dim(datass)[2])]
result=NULL
for( i in 1:length(VarNames)){
  temp=UniCox(datass,VarNames[i])
  result=rbind(result,temp)
}
result$Characteristics[result$p.value<0.05]
UniVar.sig=result[which(result$p.value<0.05),]
dir.create("CpG")
write.table(result,"CpG/single_cox_result.txt",sep="\t",quote=FALSE,row.names = F)
write.table(UniVar.sig,"CpG/single_cox_sig.txt",sep="\t",quote=FALSE,row.names = F)


#################
##-----------七、TF相关的甲基化位点的网络---------
##--获得同一模块单因素显著的转录因子和单因素显著的甲基化位点的相关性------
rm(list=ls())
library(psych)#psych包用于计算相关性、p值等信息
library(pheatmap)
library(reshape2)
##---读取基因表达和甲基化表达---#
important_path="D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count/find_important"
TF_path="D:/project/07.STAD_subtype/analysis0507_v7/04.OS/TF"
meth_sig_path="D:/project/07.STAD_subtype/analysis0507_v7/04.OS/CpG"
meth_sig=read.table(file.path(meth_sig_path,"single_cox_sig.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
TF=read.table(file.path(TF_path,"single_cox_sig.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
TF_sig=unique(TF$Characteristics)
cor=readRDS(file.path(important_path,"module_cg_gene_cor.rds"))
TF=read.table(file.path(important_path,"important_TF.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
TF=TF[TF$TF %in% TF_sig,]
module=unique(TF$module)
sig=NULL
for(i in 1:length(module)){
  temp=cor[[module[i]]]
  cor_temp=temp$sig
  cor_temp=cor_temp[cor_temp$TF %in% TF$TF,]
  cor_temp=cor_temp[cor_temp$cg %in% meth_sig$Characteristics,]
  sig=rbind(sig,cor_temp)
}
dir.create("D:/project/07.STAD_subtype/analysis0507_v7/05.cox_TF_cox_CpG_cor")
setwd("D:/project/07.STAD_subtype/analysis0507_v7/05.cox_TF_cox_CpG_cor")
write.table(sig,"TF_CpG_cor_sig.xls",sep="\t",quote=FALSE,row.names = F)
################
##--TF_regulon------
rm(list=ls())
library(dplyr)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/05.cox_TF_cox_CpG_cor")
sig=read.table("TF_CpG_cor_sig.xls",sep="\t",header=T,stringsAsFactors = F,quote="")
#sig=sig[sig$cor< -0.4,]
#TF=unique(sig$TF)
TF_cor_temp=sig %>% group_by(TF) %>% summarize(mean_cor=round(mean(cor),3),
                                               max_cor=round(max(cor),3),min_cor=round(min(cor),3),
                                               median_cor=round(median(cor),3),number=length(unique(cg)))
TF_cor_temp$regulon=paste0(TF_cor_temp$TF,"(",TF_cor_temp$number,")")
write.table(TF_cor_temp,"TF_regulon_cor.xls",sep="\t",quote=FALSE,row.names = F)
#####################
##--TF单因素cox显著的相关的甲基化位点网络的node-----
##--TF CpG network node------
rm(list=ls())
library(dplyr)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/05.cox_TF_cox_CpG_cor")
sig=read.table("TF_CpG_cor_sig.xls",sep="\t",header=T,stringsAsFactors = F,quote="")
length(unique(sig$cg))
aa=as.data.frame(table(sig$TF))
sig$abs_cor=abs(sig$cor)
write.table(sig,"TF_cox_sig_cor.txt",sep="\t",quote=FALSE,row.names = F)
TF=unique(sig$TF)
CpG=unique(sig$cg)
node=data.frame(node=c(TF,CpG),
                type=c(rep("TF",length(TF)),
                       rep("CpG",length(CpG))))
write.table(node,"TF_cox_sig_cor_node.txt",sep="\t",quote=FALSE,row.names = F)
##--TF_regulon network node-----
rm(list=ls())
library(dplyr)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/05.cox_TF_cox_CpG_cor")
TF_cor=read.table("TF_regulon_cor.xls",sep="\t",header=T,stringsAsFactors = F,quote="")
TF_cor$abs=abs(TF_cor$median_cor)
dir.create("network_regulon")
setwd("network_regulon")
write.table(TF_cor,"TF_sig_cox_regulon_cor.txt",sep="\t",quote=FALSE,row.names = F)
TF=TF_cor$TF
CpG=TF_cor$regulon
node=data.frame(node=c(TF,CpG),
                type=c(rep("TF",length(TF)),
                       rep("CpG",length(CpG))),
                size=c(TF_cor$number,rep(1,length(CpG))))
write.table(node,"TF_sig_cox_regulon_cor_node.txt",sep="\t",quote=FALSE,row.names = F)

############################
##--八、GSVA TF_regulon----
##--GSVA_score*cor分亚型------
##--根据TF对应的甲基化位点，得到每个样本的甲基化打分 GSVA-----
rm(list=ls())
meth_TF_path="D:/project/07.STAD_subtype/analysis0507_v7/03.WGCNA_count/find_important"
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/meth"
train=readRDS(file.path(data_path,"train_data.rds"))
test=readRDS(file.path(data_path,"test_data.rds"))
data=cbind(train,test)
gene_set_path="D:/project/07.STAD_subtype/analysis0507_v7/05.cox_TF_cox_CpG_cor"
##读取背景基因集合
gene_set=read.table(file.path(gene_set_path,"TF_CpG_cor_sig.xls"),header=T,sep="\t",stringsAsFactors = F)
list=split(as.matrix(gene_set)[,1], gene_set[,2])##存储的是每个免疫细胞对应的基因,构建背景基因集合
library(GSVA)
gsva_matrix=gsva(as.matrix(data),list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)##进行gsva分析
dir.create("D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation")
setwd("D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation")
write.table(gsva_matrix,"GSVA_methylation.txt",sep="\t",quote=FALSE,row.names = TRUE)##
################
##--GSVA_score*cor得分-----------
rm(list=ls())
cor_path="D:/project/07.STAD_subtype/analysis0507_v7/05.cox_TF_cox_CpG_cor"
setwd("D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation")
##---甲基化分数--##
meth_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/meth"
sample_type_meth=read.table(file.path(meth_path,"sample_type.txt"),sep="\t",quote="",header=T,stringsAsFactors = F)
sample_type_meth$sample2=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", sample_type_meth$sample)
GSVA_score=read.table("GSVA_methylation.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
GSVA_score=GSVA_score[,sample_type_meth$sample[sample_type_meth$type=="Tumor"]]
#colnames(GSVA_score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(GSVA_score))
##--TF_cor--#
TF_cor=read.table(file.path(cor_path,"TF_regulon_cor.xls"),header=T,sep="\t",stringsAsFactors = F)
TF_cg_score=NULL
for(i in 1:dim(TF_cor)[1]){
  score=GSVA_score[TF_cor$TF[i],]
  cor=TF_cor$median_cor[i]
  final_score=cor*score
  TF_cg_score=rbind(TF_cg_score,final_score)
}
write.table(TF_cg_score,"TF_cg_score.xls",sep="\t",quote=FALSE,row.names = T)
train_sample=read.table(file.path(meth_path,"train_sample.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
train_sample=train_sample[train_sample$type=="Tumor",]
test_sample=read.table(file.path(meth_path,"test_sample.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
test_sample=test_sample[test_sample$type=="Tumor",]
##
train_TF_cg_score=TF_cg_score[,train_sample$sample]
test_TF_cg_score=TF_cg_score[,test_sample$sample]
write.table(train_TF_cg_score,"train_TF_cg_score.xls",sep="\t",quote=FALSE,row.names = T)
write.table(test_TF_cg_score,"test_TF_cg_score.xls",sep="\t",quote=FALSE,row.names = T)

###--------ConsensusClusterPlus一致性聚类分亚型------#
###--------ConsensusClusterPlus一致性聚类分亚型------#
##---TF也是单因素cox显著的---------#
#--------train ConsensusClusterPlus一致性聚类分亚型------------
rm(list=ls())
data_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation"
data=read.table(file.path(data_path,"train_TF_cg_score.xls"),header=T,sep="\t",stringsAsFactors = F,row.names = 1,check.names = F)##训练集分数
library(ConsensusClusterPlus)
data=as.matrix(data) ##把数据转成R包输入格式
dir.create("ConsensusClusterPlus10-100-km-euclidean")
setwd("ConsensusClusterPlus10-100-km-euclidean")
workDir="ConsensusClusterPlus10-100-km-euclidean" ##建立放置结果的文件夹
###主语句
maxK=10
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=100,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="pdf")
results[[2]][["consensusClass"]][1:20]####查看k为2的样本分类的结果
for(i in 2:10){
  write.table(results[[i]][["consensusClass"]],paste("consensusClass",i,".txt",sep=""),quote=F,row.names=T,col.names=F,sep="\t")
} ###循环写出k从2到10的分类结果
#######################
###-------- test ConsensusClusterPlus一致性聚类分亚型------
rm(list=ls())
data_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/"
data=read.table(file.path(data_path,"test_TF_cg_score.xls"),header=T,sep="\t",stringsAsFactors = F,row.names = 1,check.names = F)##甲基化数据
library(ConsensusClusterPlus)
data=as.matrix(data) ##把数据转成R包输入格式
setwd("D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/")
dir.create("test_ConsensusClusterPlus10-100-km-euclidean")
workDir="test_ConsensusClusterPlus10-100-km-euclidean" ##建立放置结果的文件夹
setwd("test_ConsensusClusterPlus10-100-km-euclidean")
###主语句
maxK=10
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=100,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="km",
                             distance="euclidean",
                             seed=555,
                             plot="pdf")
###
results[[2]][["consensusClass"]][1:20]####查看k为2的样本分类的结果
for(i in 2:10){
  write.table(results[[i]][["consensusClass"]],paste("consensusClass",i,".txt",sep=""),quote=F,row.names=T,col.names=F,sep="\t")
} ###循环写出k从2到10的分类结果

###

##---将训练集的临床结果和分类合并在一起分两类------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/ConsensusClusterPlus10-100-km-euclidean")
class=read.table("consensusClass2.txt",sep="\t",header=F,stringsAsFactors = F,quote="")
colnames(class)=c("sample","class")
class$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", class$sample)
clinical_path="D:/project/07.zhouY_project/01.analysis0422/00.data/clinical/"
clinical=read.table(file.path(clinical_path,"TCGA-STAD-tumor.survival.tsv"),sep="\t",header=T,stringsAsFactors = F,quote="")
survival=merge(clinical,class,by="sample")
write.table(survival,"survival.txt",sep="\t",quote=FALSE,row.names = F)
##---绘制训练集不同亚型的生存曲线 分两类-----
rm(list=ls())
library(survival)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/ConsensusClusterPlus10-100-km-euclidean")
survival1=read.table("survival.txt",header=T,sep="\t",stringsAsFactors = F)
surv=survfit(Surv(survival1$OS.time, survival1$OS) ~ survival1$class, data = survival1) ##填入的数据分别是生存状态，时间，类别
##---使用survminer绘制生存曲线--#
library(survminer)
pdf("survminer.pdf",width=7,height=6)
print(ggsurvplot(surv,
                 #legend.labs=c("cluster1", "cluster2"), #图注
                 legend= c(0.9,0.5),#图例位置
                 #conf.int = TRUE,
                 xlab = "Days",
                 pval=TRUE,#添加p值
                 pval.method = TRUE,
                 #pval.method.size = 2.5,
                 legend.title = "Cluster",
                 group.table=T,
                 cumevents=F,
                 risk.table = T,
                 title="Survival of each set"))
dev.off()



##---将测试集临床结果和分类合并在一起 分两类-------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/test_ConsensusClusterPlus10-100-km-euclidean")
class=read.table("consensusClass2.txt",sep="\t",header=F,stringsAsFactors = F,quote="")
colnames(class)=c("sample","class")
class$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", class$sample)
clinical_path="D:/project/07.zhouY_project/01.analysis0422/00.data/clinical/"
clinical=read.table(file.path(clinical_path,"TCGA-STAD-tumor.survival.tsv"),sep="\t",header=T,stringsAsFactors = F,quote="")
survival=merge(clinical,class,by="sample")
write.table(survival,"test_survival.txt",sep="\t",quote=FALSE,row.names = F)
##---绘制测试集不同亚型的生存曲线分两类-----
rm(list=ls())
library(survival)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/test_ConsensusClusterPlus10-100-km-euclidean")
survival1=read.table("test_survival.txt",header=T,sep="\t",stringsAsFactors = F)
surv=survfit(Surv(survival1$OS.time, survival1$OS) ~ survival1$class, data = survival1) ##填入的数据分别是生存状态，时间，类别
##---使用survminer绘制生存曲线--#
library(survminer)
pdf("test_survminer.pdf",width=7,height=6)
print(ggsurvplot(surv,
                 #legend.labs=c("cluster1", "cluster2"), #图注
                 legend= c(0.9,0.5),#图例位置
                 #conf.int = TRUE,
                 xlab = "Days",
                 pval=TRUE,#添加p值
                 pval.method = TRUE,
                 #pval.method.size = 2.5,
                 legend.title = "Cluster",
                 group.table=T,
                 cumevents=F,
                 risk.table = T,
                 title="Survival of each set"))
dev.off()

##---rename 2类的----
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/")
train=read.table("ConsensusClusterPlus10-100-km-euclidean/consensusClass2.txt",sep="\t",header = F,stringsAsFactors = F,quote="")
colnames(train)=c("sample","cluster")
test=read.table("test_ConsensusClusterPlus10-100-km-euclidean/consensusClass2.txt",sep="\t",header = F,stringsAsFactors = F,quote="")
colnames(test)=c("sample","cluster")
train$cluster[train$cluster==2]="cluster1"
train$cluster[train$cluster==1]="cluster2"
test$cluster[test$cluster==1]="cluster1"
test$cluster[test$cluster==2]="cluster2"
dir.create("cluster_divide2")
setwd("cluster_divide2")
write.table(train,"train_cluster.xls",sep="\t",quote=FALSE,row.names = F)
write.table(test,"test_cluster.xls",sep="\t",quote=FALSE,row.names = F)
##---survival 2类的----
rm(list=ls())
library(survival)
library(survminer)
# ##---survival function---##
# surv_add_CI <- function(data1,title1){
#   data1=data1
#   title1=title1
#   fits_OS <- survfit(Surv(OS.time,OS) ~ class, data = data1)
#   mm<-survival:::survmean(fits_OS,rmean = "none")
#   res_cox<-coxph(Surv(OS.time,OS) ~ class, data = data1)
#   summary(res_cox)$conf.int 
#   HR<-round(summary(res_cox)$conf.int[1],2)   
#   CI.L<-round(summary(res_cox)$conf.int[3],2)   
#   CI.H<-round(summary(res_cox)$conf.int[4],2)
#   names(fits_OS$strata)<-gsub("class=","",names(fits_OS$strata))
#   pvalue <- pairwise_survdiff(Surv(OS.time,OS) ~ class, data = data1)
#   all_p2 <- round(signif(pvalue$p.value,3),3)
#   p <- ggsurvplot(fits_OS,
#                   legend.title = element_blank(),
#                   legend=c(0.8,0.9), 
#                   #break.time.by=5,
#                   title=paste0("Overall survival of ",title1), 
#                   #pval = T,# 添加P值
#                   #pval.method=T,
#                   size = 0.6,
#                   risk.table = TRUE,# 添加风险表
#                   #risk.table.height=0.30,
#                   tables.height = 0.2,
#                   #tables.theme = theme_cleantable(),
#                   surv.median.line = "hv",# 增加中位生存时间
#                   #conf.int = T,# 增加置信区间
#                   #add.all = TRUE, # 添加总患者生存曲线
#                   censor.size=5,
#                   palette = c("#E7B800", "#2E9FDF"),
#                   ggtheme=theme_survminer(font.main = c(12,"plain","black"),
#                                           font.x = c(12,"plain","black"),
#                                           font.y = c(12,"plain","black"),
#                                           font.caption = c(12,"plain","black"),
#                                           font.tickslab = c(12,"plain","black"))) 
#   p$plot <- p$plot +xlab("Time(days)")+ ylab("Survival Probability")+ 
#     annotate("text",x=25,y=0.5,label=paste0("Log-rank: p=",all_p2,"\nHR = ",HR," (95% CI: ",CI.L,"-",CI.H,")"),hjust=0,size=5) +
#     annotate("text",x=25,y=0.72,label=paste0(names(fits_OS$strata)[1],"\n(median ",mm$matrix[,"median"][1]," months, 95% CI: ",
#                                              mm$matrix[,"0.95LCL"][1],"-",mm$matrix[,"0.95UCL"][1],")"),hjust=0,size=5)+
#     annotate("text",x=25,y=0.6,label=paste0(names(fits_OS$strata)[2],"\n(median ",mm$matrix[,"median"][2]," months, 95% CI: ",
#                                             mm$matrix[,"0.95LCL"][2],"-",mm$matrix[,"0.95UCL"][2],")"),hjust=0,size=5)+
#     theme(plot.title = element_text(hjust = 0.5),
#           legend.text = element_text(size=9),
#           legend.background = element_blank())
#   return(p)
# }

#---train survival----#
setwd("D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/")
train=read.table("ConsensusClusterPlus10-100-km-euclidean/survival.txt",header=T,sep="\t",stringsAsFactors = F)
train$class[train$class==2]="cluster1"
train$class[train$class==1]="cluster2"
title_train="Train Set"
##---test survival---##
test=read.table("test_ConsensusClusterPlus10-100-km-euclidean/test_survival.txt",header=T,sep="\t",stringsAsFactors = F)
test$class[test$class==1]="cluster1"
test$class[test$class==2]="cluster2"
title_test="Test Set"
p_list=list()
for(i in 1:2){
  if(i==1){
    data1=train
    title1=title_train
  }else{
    data1=test
    title1=title_test
  }
  fits_OS <- survfit(Surv(OS.time,OS) ~ class, data = data1)
  mm<-survival:::survmean(fits_OS,rmean = "none")
  res_cox<-coxph(Surv(OS.time,OS) ~ class, data = data1)
  summary(res_cox)$conf.int 
  HR<-round(summary(res_cox)$conf.int[1],2)   
  CI.L<-round(summary(res_cox)$conf.int[3],2)   
  CI.H<-round(summary(res_cox)$conf.int[4],2)
  names(fits_OS$strata)<-gsub("class=","",names(fits_OS$strata))
  pvalue <- pairwise_survdiff(Surv(OS.time,OS) ~ class, data = data1)
  all_p2 <- round(signif(pvalue$p.value,3),3)
  p <- ggsurvplot(fits_OS,
                  legend.title = element_blank(),
                  legend=c(0.8,0.9), 
                  #break.time.by=5,
                  title=paste0("Overall survival of ",title1), 
                  #pval = T,# 添加P值
                  #pval.method=T,
                  size = 0.6,
                  risk.table = TRUE,# 添加风险表
                  #risk.table.height=0.30,
                  tables.height = 0.2,
                  #tables.theme = theme_cleantable(),
                  surv.median.line = "hv",# 增加中位生存时间
                  #conf.int = T,# 增加置信区间
                  #add.all = TRUE, # 添加总患者生存曲线
                  censor.size=5,
                  palette = c("#61b3de","#ffee93"),
                  ggtheme=theme_survminer(font.main = c(12,"plain","black"),
                                          font.x = c(12,"plain","black"),
                                          font.y = c(12,"plain","black"),
                                          font.caption = c(12,"plain","black"),
                                          font.tickslab = c(12,"plain","black"))) 
  p$plot <- p$plot +xlab("Time(days)")+ ylab("Survival Probability")+ 
    annotate("text",x=25,y=0.05,label=paste0("Log-rank: p=",all_p2,"\nHR = ",HR," (95% CI: ",CI.L,"-",CI.H,")"),hjust=0,size=3.5) +
    annotate("text",x=25,y=0.35,label=paste0(names(fits_OS$strata)[1],"\n(median ",mm$matrix[,"median"][1]," days, 95% CI: ",
                                             mm$matrix[,"0.95LCL"][1],"-",mm$matrix[,"0.95UCL"][1],")"),hjust=0,size=3.5)+
    annotate("text",x=25,y=0.2,label=paste0(names(fits_OS$strata)[2],"\n(median ",mm$matrix[,"median"][2]," days, 95% CI: ",
                                            mm$matrix[,"0.95LCL"][2],"-",mm$matrix[,"0.95UCL"][2],")"),hjust=0,size=3.5)+
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size=9),
          legend.background = element_blank())
  p_list[[i]]=p
}
library(survminer)
setwd("cluster_divide2")
pdf("surv_add_CI.pdf",width=12,height=6)
arrange_ggsurvplots(p_list,print=TRUE,ncol=2,risk.table.height=.25)
dev.off()
pdf("surv_add_CI_train.pdf",width=6,height=6)
p_list[[1]]
dev.off()
pdf("surv_add_CI_test.pdf",width=6,height=6)
p_list[[2]]
dev.off()
##---PCA 两类-------
##--train PCA 两类---------
rm(list=ls())
data_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/"
data=read.table(file.path(data_path,"train_TF_cg_score.xls"),header=T,sep="\t",stringsAsFactors = F,row.names = 1,check.names = F)##
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
cluster=read.table(file.path(cluster_path,"train_cluster.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
rownames(cluster)=cluster$sample
sample=intersect(cluster$sample,colnames(data))
data=data[,sample]
cluster=cluster[sample,]
data.pca<-princomp(t(data),cor=T,scores=T) 
library(ggbiplot)
setwd(cluster_path)
pdf("train_PCA.pdf",width = 6,height=6)
ggbiplot(data.pca, obs.scale = 1, var.scale = 1,
         groups = cluster$cluster, ellipse = FALSE, var.axes = F)+
  scale_color_manual(name="Cluster", values=c("#61b3de","#ffee93"))

dev.off()
#install.packages("devtools", repo="http://cran.us.r-project.org")
#library(devtools)
#install_github("vqv/ggbiplot")
#library(remotes)
#remotes::install_git("https://gitee.com/git_lzl/ggbiplot.git",force = TRUE)
##-----test PCA 两类-----
rm(list=ls())
data_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/"
data=read.table(file.path(data_path,"test_TF_cg_score.xls"),header=T,sep="\t",stringsAsFactors = F,row.names = 1,check.names = F)##甲基化数据
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
cluster=read.table(file.path(cluster_path,"test_cluster.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
rownames(cluster)=cluster$sample
sample=intersect(cluster$sample,colnames(data))
data=data[,sample]
cluster=cluster[sample,]
data.pca<-princomp(t(data),cor=T,scores=T) 
library(ggbiplot)
setwd(cluster_path)
pdf("test_PCA.pdf",width = 6,height=6)
ggbiplot(data.pca, obs.scale = 1, var.scale = 1,
         groups = cluster$cluster, ellipse = FALSE, var.axes = F)+
  scale_color_manual(name="Cluster", values=c("#61b3de","#ffee93"))

dev.off()

##---ggplot2 画PCA 两类--------
rm(list=ls())
data_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/"
train_score=read.table(file.path(data_path,"train_TF_cg_score.xls"),header=T,sep="\t",stringsAsFactors = F,row.names = 1,check.names = F)##
test_score=read.table(file.path(data_path,"test_TF_cg_score.xls"),header=T,sep="\t",stringsAsFactors = F,row.names = 1,check.names = F)##
cluster_path="D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2"
train_cluster=read.table(file.path(cluster_path,"train_cluster.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
rownames(train_cluster)=train_cluster$sample
test_cluster=read.table(file.path(cluster_path,"test_cluster.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
rownames(test_cluster)=test_cluster$sample
PCA_ggplot <- function(cluster,data,title,size){
  cluster$cluster[cluster$cluster=="immune_enriched"]="cluster1"
  cluster$cluster[cluster$cluster=="immune_poor"]="cluster2"
  sample=intersect(cluster$sample,colnames(data))
  data=data[,sample]
  cluster=cluster[sample,]
  #主成分计算
  pca_data=prcomp(t(data), scale = TRUE)
  #查看合适主成分个数
  screeplot(pca_data, type = "lines")
  # summary(pca_data)
  # #查看行名，确认是否为36个样品的名称
  # rownames(pca_data$x)
  #提取PC1的百分比
  x_per=round(summary(pca_data)$importance[2, 1]*100, 1)
  #提取PC2的百分比
  y_per=round(summary(pca_data)$importance[2, 2]*100, 1)
  #按照样品名称添加组并且合并
  df_sample=data.frame(sample=rownames(pca_data$x), pca_data$x) %>%
    left_join(cluster, by = "sample")
  #绘图
  plot_1 <- ggplot(df_sample, aes(x = PC1, y = PC2)) +
    geom_point(aes(color=cluster),size=size) +
    xlab(paste("PC1","(", x_per,"%)",sep=" ")) +
    ylab(paste("PC2","(", y_per,"%)",sep=" ")) +
    labs(title=title)+
    scale_color_manual(name="Cluster", values=c("#61b3de","#ffee93"))+
    theme_bw()
  return(plot_1)
  
}
plot_train <- PCA_ggplot(train_cluster,train_score,"Train Datesets",size=2)
plot_test <- PCA_ggplot(test_cluster,test_score,"Test Datesets",size=2.5)
plot2 <- list(plot_train,plot_test)
library(cowplot)
pdf("train_PCA_v2.pdf",width=5,height=5)
plot_train
dev.off()
pdf("test_PCA_v2.pdf",width=5,height=5)
plot_test
dev.off()
pdf("PCA_v2.pdf",width=8,height=3.5)
plot_grid(plotlist=plot2,ncol=2,rel_widths = c(1.02,1)
)
dev.off()
#############

##---TF regulon pheatmap train test------
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
train=read.table(file.path(cluster_path,"train_cluster.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
test=read.table(file.path(cluster_path,"test_cluster.xls"),sep="\t",header=T,stringsAsFactors = F,quote="")
train$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", train$sample)
rownames(train)=train$sample
train$cluster=factor(train$cluster,levels=c("cluster1","cluster2"))
test$sample=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", test$sample)
rownames(test)=test$sample
test$cluster=factor(test$cluster,levels=c("cluster1","cluster2"))
train_cluster=train
test_cluster=test
network_path="D:/project/07.STAD_subtype/analysis0507_v7/05.cox_TF_cox_CpG_cor/network_regulon"
network=read.table(file.path(network_path,"TF_sig_cox_regulon_cor.txt"),sep='\t',header=T,stringsAsFactors = F,quote="")
network=network[order(network$number,network$abs,decreasing = T),]
network1=network[,c("TF","regulon")]
pheatmap_clincial <- function(score,TF,cluster,network1){
  data=score[TF$Characteristics,]
  cluster=cluster[order(cluster$cluster),]
  data=data[,cluster$sample]
  data=merge(network1,data,by.x="TF",by.y="row.names")
  rownames(data)=data$regulon
  data=data[,-c(1,2)]
  data=data[network1$regulon,]
  ##---color---##
  clustercolor <- c( "#61b3de","#ffee93") 
  names(clustercolor) <- c("cluster1","cluster2") #类型颜色
  ann_colors=list(cluster=clustercolor)
  annotation_col=data.frame(cluster=cluster$cluster)
  row.names(annotation_col)=cluster$sample
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
p_list=pheatmap_clincial(train_score,TF,train_cluster,network1)
p_list_test=pheatmap_clincial(test_score,TF,test_cluster,network1)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/06.GSVA_methylation/cluster_divide2")
pdf("train_pheatmap.pdf",width=8,height=6.5)
p_list[[1]]
dev.off()
pdf("train_pheatmap_no_cluster.pdf",width=8,height=6.5)
p_list[[2]]
dev.off()
pdf("test_pheatmap.pdf",width=8,height=6.5)
p_list_test[[1]]
dev.off()
pdf("test_pheatmap_no_cluster.pdf",width=8,height=6.5)
p_list_test[[2]]
dev.off()
#################