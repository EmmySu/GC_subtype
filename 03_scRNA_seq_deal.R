##----九、利用单细胞定义不同亚型------
##----合并GSE167297的肿瘤样本----
rm(list=ls())
library(Seurat)
library(cowplot)
library(dplyr)
library(glmGamPoi)
scRNA_path=("D:/project/07.STAD_subtype/analysis0507_v7/00.data/GEO_scRNA-seq/GSE167297/")
file=list.files(scRNA_path)
file=file[-grep("Normal",file)]
get_MT_GEO <- function(file){
  infile <- paste(scRNA_path,file,sep="")
  samplename <- str_split_fixed(file,"_",2)[,1]
  Data <- read.table(infile, sep="\t", header=T, row.names=1) 
  data <- CreateSeuratObject(counts = Data, project = samplename, min.cells = 3, min.features = 200)
  ###----2.细胞质控，数据过滤---##
  ## 计算线粒体占比和rRNA
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")# 计算线粒体占比,MT-开头的基因认为是线粒体基因，鼠源的需要替换为mt
  #data[["percent.pHB"]] <- PercentageFeatureSet(data, pattern = "^HBA|^HBB")##血红蛋白
  #data[["percent.pRP"]] <- PercentageFeatureSet(data, pattern = "^RPS|^RPL")##核糖体基因
  #nFeature_RNA代表每个细胞测到的基因数目，nCount代表每个细胞测到所有基因的表达量之和，percent.mt代表测到的线粒体基因的比
  p=VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"
                               # "percent.pHB",
                               #"percent.pRP"
  ), ncol = 3)
  sample_MT <- list(data=data,p=p)
  return(sample_MT)
}

get_object_GEO <- function(data,samplename,nFeature_low,nFeature_high,mt){
  # ###--这里是要过滤掉多重细胞--##
  # scrubletdir=paste("D:/linux/project/01.20220825STAD/multipletsCell/",samplename,sep="")
  # multipleCell <- read.table(paste(scrubletdir,"/", samplename, ".predicted_doublets_barcode.txt", sep = ''), header = T, sep = "\t", check.names = F)
  # multipleCell <- subset(multipleCell, Double == "False")
  # data <- subset(data, cells = multipleCell$Barcode)
  ####--过滤细胞--#
  data <- subset(data, subset = nFeature_RNA > nFeature_low & nFeature_RNA < nFeature_high & 
                   percent.mt < mt 
                 #& percent.pRP<pRP
  ) #根据前面的图来确定标准
  threshold <- data.frame(sample=samplename,nFeature_low=nFeature_low,nFeature_high=nFeature_high,mt=mt
                          #,pRP=pRP
  )
  object <- list(data=data,threshold=threshold)
  return(object)
}
threshold=NULL
datalist=list()
sample=NULL
setwd("D:/project/07.STAD_subtype/analysis0507_v7/07.scRNA-seq")
for(i in  1:length(file)){
  data_MT <- get_MT_GEO(file[i])
  data <- data_MT$data
  p <- data_MT$p
  sample_GEO <- str_split_fixed(file[i],"_",2)[,1]
  sampledir <- sample_GEO
  if (! dir.exists(sampledir)){
    dir.create(sampledir)
  }
  pdf(paste(sampledir,"/QCmetrics.pdf",sep=""),onefile = F)
  print(p)
  dev.off()
  
  data <- get_object_GEO(data,sample_GEO,500,4000,10)
  threshold_temp <- data.frame(sample=sample_GEO,nFeature_low=500,nFeature_high=4000,mt=10)
  threshold=rbind(threshold,threshold_temp)
  sample=c(sample,sample_GEO)
  data=data$data
  #data@meta.data$orig.ident
  datalist[[i]]=data
}
names(datalist)=sample
print(datalist)
immune.combined <- merge(datalist[[1]],datalist[2:length(datalist)], add.cell.ids =  sample)
# FeatureScatter(immune.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# VlnPlot(immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/07.scRNA-seq")
dir.create("merge")
saveRDS(immune.combined,"merge/Seurat_just_merge.rds")
write.table(threshold,"sample_threshold.xls",sep="\t",quote=FALSE,row.names = F)
##############################
###---harmony去批次-----------
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/07.scRNA-seq")
library(Seurat)
library(ggsci)
library(scales)
library(ggplot2)
library(patchwork)
data=readRDS("merge/Seurat_just_merge.rds")
unique(data@meta.data$orig.ident)
###---QC plot----#
sample_col=c("#BC3C29FF","#E18727FF","#20854EFF", "#7876B1FF" ,"#6F99ADFF","#FB8072", "#EDC948","#B2DF8A",  "#BEBADA","#76B7B2")
theme.set2=theme(axis.title.x=element_blank())
plot.features=c("nFeature_RNA","nCount_RNA","percent.mt")
group="orig.ident"
plots=list()
plots2=list()
for(i in seq_along(plot.features)){
  plots[[i]]=VlnPlot(data,group.by = group,
                     #pt.size=0,
                     features=plot.features[i],
                     cols=sample_col)+theme.set2+NoLegend()
  plots2[[i]]=VlnPlot(data,group.by = group,pt.size=0,
                      features=plot.features[i],
                      cols=sample_col)+theme.set2+NoLegend()
  
}
violin=wrap_plots(plots=plots,nrow=2)
violin2=wrap_plots(plots=plots2,nrow=2)
ggsave("merge/vlnplot_merge_add_points.pdf",plot=violin,width=10,height=8)
ggsave("merge/vlnplot_merge.pdf",plot=violin2,width=10,height=6)
##--relationship plot---#
plot_cor=list()
plot_cor[[1]]=FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt",cols=sample_col)+NoLegend()
point=wrap_plots(plots=plot_cor,nrow=1)
ggsave("merge/relationship_merge.pdf",plot=point,width=5,height=5)
##---NormalizeData--#
data=SCTransform(data, method = "glmGamPoi", vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = FALSE)
data=FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10=head(VariableFeatures(data), 10)
# plot variable features with and without labels
plot1=VariableFeaturePlot(data)
plot2=LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave("merge/huoshantu.pdf",plot=plot2,width=7.5,height=5)
##---cell cycle---#
g2m_genes=cc.genes$g2m.genes
g2m_genes=CaseMatch(search = g2m_genes, match = rownames(data))
s_genes=cc.genes$s.genes
s_genes=CaseMatch(search = s_genes, match = rownames(data))
data=CellCycleScoring(data, g2m.features = g2m_genes, s.features = s_genes, set.ident = FALSE)
data=SCTransform(data, method = "glmGamPoi", vars.to.regress = c("S.Score", "G2M.Score"))
##---PCA--#
data[["RNA"]]@counts[1:4,1:4]
data=RunPCA(data, features = VariableFeatures(data), verbose = FALSE)
p=DimPlot(data, reduction = "pca")
ggsave("merge/PCA.pdf",plot=p,width=7.5,height=5)
# 应该选多少个主成分进行后续分析
pca=ElbowPlot(data,ndims=50)
ggsave("merge/PCA_choose.pdf",plot=pca,width=7.5,height=5)
saveRDS(data,"merge/Seurat_before_harmony.rds")
##---harmony---#
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/07.scRNA-seq")
library(harmony)
data=readRDS("merge/Seurat_before_harmony.rds")
data=RunHarmony(data,group.by="orig.ident",assay.use="SCT",max.inter.harmony=20)
## max.inter.harmony是设置迭代次数，默认是10
##RunHarmony函数中有个lambda参数，默认值是1，决定了harmony的整合的力度，lambda值调小，整合力度变大，调整范围一般是0.5-2之间
saveRDS(data,"merge/Seurat_after_harmony.rds")
##----分亚群-----------------
##----cluster-----------
###--harmony---#
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/07.scRNA-seq")
library(Seurat)
library(ggthemes)
library(dplyr)
library(ggplot2)
data=readRDS("merge/Seurat_after_harmony.rds")
pca=ElbowPlot(data,ndims=50)
pc.num=1:30
data=RunTSNE(data, reduction="harmony", dims=pc.num) %>% 
  RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num) %>% 
  FindClusters(resolution=0.5)
length(levels(Idents(data)))
table(data@active.ident) 
dir.create("harmony_0.5_30")
setwd("harmony_0.5_30")
sample_col=c("#BC3C29FF","#E18727FF","#20854EFF", "#7876B1FF" ,"#6F99ADFF","#FB8072", "#EDC948","#B2DF8A",  "#BEBADA","#76B7B2")
##color
p=DimPlot(data,reduction = "umap",label=T)
p1=DimPlot(data,reduction = "umap",group.by="orig.ident",split.by="orig.ident",ncol=3,cols=sample_col)+NoLegend()
p11=DimPlot(data,reduction = "umap",group.by="orig.ident",cols=sample_col)
ggsave("UMAP_Samples_harmony_no_split.pdf",p11,width=7.5,height=5)
ggsave("UMAP_Samples_harmony.pdf",p1,width=9,height=9)
ggsave("UMAP_cluster_harmony.pdf",p,width=7.5,height=5)
p2=DimPlot(data,reduction = "tsne",group.by="orig.ident",split.by="orig.ident",ncol=3,cols=sample_col)+NoLegend()
p22=DimPlot(data,reduction = "tsne",group.by="orig.ident",cols=sample_col)
p3=DimPlot(data,reduction = "tsne",label=T)
ggsave("TSNE_Samples_harmony_no_split.pdf",p22,width=7.5,height=5)
ggsave("TSNE_Samples_harmony.pdf",p2,width=9,height=9)
ggsave("TSNE_cluster_harmony.pdf",p3,width=7.5,height=5)
saveRDS(data,"Seurat.rds")
##---cluster marker---#
clus=data.frame(data@active.ident)###cluster 的结果
colnames(clus)=c("Cluster")
dir.create("Markers")
write.csv(clus,'Markers/Cluster.csv',quote=F)
##---find markers---#
data.markers=FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) ##里面还有一个参数test.use 可以选取不同的方法 有wilcox,MAST等，默认是wilcox
write.table(data.markers,'Markers/Markers.xls',sep='\t',quote=F,col.names=T, row.names=F)
Htop10=data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
print(Htop10$gene)
Htop5=data.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) 
Htop2=data.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) 
pdf('Markers/Heatmap_top2.pdf',onefile=F,width = 10,height=8)
DoHeatmap(data, features = Htop2$gene) + NoLegend()
dev.off()
pdf('Markers/Heatmap_top5.pdf',onefile=F)
DoHeatmap(data, features = Htop5$gene) + NoLegend()
dev.off()
###########
########----------singleR 注释---------
##---使用singleR HCL human 注释---#
rm(list=ls())
setwd("D:/project/07.STAD_subtype/analysis0507_v7/07.scRNA-seq/harmony_0.5_30")
Sample="STAD"
Database="HP"
type="umap"
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleR)
library(tibble)
options(future.globals.maxSize= 1300 * 1024^2)
### read 10X data
#---- five human ref databse ----##
if (Database == "HP"){
  load("E:/00.database/SingleR/HumanPrimaryCellAtlas_hpca.se_human.RData")
  REF=hpca.se
}
if (Database == "BL"){
  load("E:/00.database/database/SingleR/BlueprintEncode_bpe.se_human.RData")
  REF=bpe.se
}
if (Database == "DICE"){
  load("E:/00.database/SingleR/DatabaseImmuneCellExpression_immu_human.RData")
  REF=immu
}
if (Database == "NH"){
  load("E:/00.database/SingleR/NovershternHematopoietic_nh_human.RData")
  REF=nh
}
if (Database == "MI"){
  load("E:/00.database/SingleR/MonacoImmune_mi_human.RData")
  REF=mi
}
#---- two mouse ref databse ---##
if (Database == "IG_M"){
  load("E:/00.database/SingleR/ImmGenData_immu.se_mouse.RData")
  REF=immu.se
}
if (Database == "MR_M"){
  load("E:/00.database/SingleR/MouseRNAseq_mr_mouse.RData")
  REF=mr
}
if (Database != "HP" & Database != "BL" & Database != "DICE" & Database != "NH" & Database != "MI" & Database != "IG_M" & Database !="MR_M" & Database !="HCL"){
  file=paste("E:/00.database/HCL/",Database,".RData",sep="")
  print(file)
  load(file)
  REF=ref	
}
ctrl=readRDS("Seurat.rds")
ctrl@meta.data$cell.type=Idents(ctrl)
test=as.SingleCellExperiment(ctrl)
## Annot
Anno=SingleR(test = test,
             ref = REF,
             labels = REF$label.main,
             method = "cluster",
             cluster = test$cell.type)
Anno$cluster=rownames(Anno)

### extract anno 
fin=Anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
new.cluster.ids=fin$labels
names(new.cluster.ids)=fin$cluster
ctrl=RenameIdents(ctrl, new.cluster.ids)
dir.create("singleR")
setwd("singleR")
dir.create("HP")
setwd("HP")
write.table(fin,file="cluster_after.xls",sep="\t",row.names = FALSE,quote=F)
#### new umap and diff genes
if (type == "umap") {
  p1=DimPlot(ctrl, reduction = "umap", pt.size = 0.5,label=TRUE)
  pdf(paste(Sample,"_UMAP_Anno.pdf",sep=""),width=10,height=6)
  print(p1)
  # + ggtitle(label = "UMAP")
  dev.off()
  png(paste(Sample,"_UMAP_Anno.png",sep=""))
  # + ggtitle(label = "UMAP")
  print(p1)
  dev.off()
  
  p2=DimPlot(ctrl, reduction = "tsne", pt.size = 0.5,label=TRUE)
  pdf(paste(Sample, "_TSNE_Anno.pdf",sep=""),width=10,height=6)
  print(p2)
  dev.off()
}
if (type == "tsne") {
  p1=DimPlot(ctrl, reduction = "tsne", pt.size = 0.5,label=TRUE)
  pdf(paste(Sample,"_TSNE_Anno.pdf",sep=""),width=10,height=6)
  print(p1)
  # + ggtitle(label = "UMAP")
  dev.off()
  png(paste(Sample,"_TSNE_Anno.png",sep=""))
  print(p1)
  # + ggtitle(label = "UMAP")
  dev.off()
  
  p2=DimPlot(ctrl, reduction = "umap", pt.size = 0.5,label=TRUE)
  pdf(paste(Sample,"_UMAP_Anno.pdf",sep=""),width=10,height=6)
  print(p2)
  dev.off()
}

ctrl@meta.data$cell_type=Idents(ctrl)
saveRDS(ctrl, file = paste(Sample,"_Seurat_Anno_singleR.rds",sep=""))
pdf(paste(Sample,"_SingleR_heatmap.pdf",sep=""),width=8,height=5)
plotScoreHeatmap(Anno)
dev.off()
tiff(paste(Sample,"_SingleR_heatmap.tiff",sep=""))
plotScoreHeatmap(Anno)
dev.off()
write.csv(Idents(ctrl),file = paste(Sample,"_cluster_Anno.csv",sep=""),quote=F)
###------------#
rm(list=ls())
##--细胞注释------
setwd("D:/project/07.STAD_subtype/analysis0507_v7/07.scRNA-seq/harmony_0.5_30")
#remotes::install_github("lyc-1995/MySeuratWrappers")
library(MySeuratWrappers)
library(Seurat)
library(RColorBrewer)
library(ggsci)
library(ggthemes)
library(ggplot2)
library(cowplot)
data=readRDS("Seurat.rds")
metadata=data@meta.data
gene2=c("CD3D","CD2","CD3E","IL7R",#T细胞
        "GNLY","NKG7","GZMA","GZMB",#NK细胞
        "CD79A","MS4A1","MZB1",#B
        "EPCAM","KRT18",#epi
        "CD68","CD86",#myeloid
        "VWF","ENG",#endothelial
        "THY1","COL1A2",#ﬁbroblasts
        "CPA3","TPSAB1" #mast
        # "FCER1A","LYZ",#dendritric
        # "C1QA","CD14"#Mono&#Macrophage
        
        
)
p=VlnPlot(data, features = gene2,pt.size = 0,stacked = T)
pdf("vlnplot.pdf",width=5,height=6)
p
dev.off()
p2=DotPlot(data,features = gene2, cols = c("#fddbc7", "red"), dot.scale = 6) + RotatedAxis()+coord_flip()
pdf("DotPlot.pdf",width=6,height=6)
p2
dev.off()
# FeaturePlot(data,features = gene2,pt.size = 0,ncol=2,reduction = "tsne")
# FeaturePlot(data,features = gene2,pt.size = 0,ncol=2,reduction = "umap")

##---重注释--#
immune.combined.new=RenameIdents(data, `0` = "B cells", `1` = "T cells", `2` = "B cells",`7` = "Endothelial cells",
                                 `3`="T cells",`4` = "NK cells", `5` = "Myeloid cells",`6`="T cells",
                                 `8` = "Fibroblasts cells",`9`="Myeloid cells",`10`="Epithelial cells",
                                 `11` = "Mast cells", `12` = "T cells",`13`="T cells")
immune.combined.new@meta.data$cell_type_immune=Idents(immune.combined.new)
dir.create("markers_by_self")
setwd("markers_by_self")
saveRDS(immune.combined.new, file = "Seurat_Anno_cells_by_self.rds")
clus=data.frame(immune.combined.new@active.ident)
colnames(clus)=c("Cluster")
write.table(clus,'Cluster.xls', sep="\t",quote=FALSE)
new.markers=FindAllMarkers(immune.combined.new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(new.markers,'Cluster_markers.csv',quote=F,row.names = F)
###----plot---#
rm(list=ls())
library(Seurat)
library(ggplot2)
library(scales)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/07.scRNA-seq/harmony_0.5_30/markers_by_self")
data.combined.new=readRDS("Seurat_Anno_cells_by_self.rds")
palettes=ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
mycolors = c("#FB8072","#80B1D3","#B2DF8A", "#0072B5FF",   "#BEBADA",
             "#E18727FF","#20854EFF","#FFFFB3")
sample=unique(data.combined.new@meta.data$orig.ident)
sample_col=c("#BC3C29FF","#E18727FF","#20854EFF", "#7876B1FF" ,"#6F99ADFF","#FB8072", "#EDC948","#B2DF8A",  "#BEBADA","#76B7B2")
cluster=unique(data.combined.new@meta.data$cell_type_immune)
p1_after=DimPlot(data.combined.new, reduction = "umap", label = FALSE,cols=mycolors[1:length(cluster)])
p2_after=DimPlot(data.combined.new, reduction = "umap", group.by ="orig.ident",label=FALSE,cols=sample_col[1:length(sample)])
p3_after=DimPlot(data.combined.new, reduction = "tsne", label = FALSE,cols=mycolors[1:length(cluster)])
p4_after=DimPlot(data.combined.new, reduction = "tsne", group.by ="orig.ident",label=FALSE,cols=sample_col[1:length(sample)])
pdf('Cluster-UMAP.pdf',onefile=F,width=6.5,height=5)
plot_grid(p1_after)
dev.off()
pdf('Samples-UMAP.pdf',onefile=F,width=7.5,height=5)
plot_grid(p2_after)
dev.off()
pdf('Cluster-TSNE.pdf',onefile=F,width=6.5,height=5)
plot_grid(p3_after)
dev.off()
pdf('Samples-TSNE.pdf',onefile=F,width=7.5,height=5)
plot_grid(p4_after)
dev.off()


##---查看重要的TF在细胞里的表达-------
rm(list=ls())
library(MySeuratWrappers)
library(Seurat)
library(RColorBrewer)
library(ggsci)
library(ggthemes)
library(ggplot2)
library(cowplot)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/07.scRNA-seq/harmony_0.5_30/markers_by_self")
data.combined.new=readRDS("Seurat_Anno_cells_by_self.rds")
mycolors = c("#FB8072","#80B1D3","#B2DF8A", "#0072B5FF",   "#BE0BADA",
             "#E18727FF","#20854EFF","#FFFFB3")
sample=unique(data.combined.new@meta.data$orig.ident)
sample_col=c("#BC3C29FF","#E18727FF","#20854EFF", "#7876B1FF" ,"#6F99ADFF","#FB8072", "#EDC948","#B2DF8A",  "#BEBADA","#76B7B2")
cluster=unique(data.combined.new@meta.data$cell_type_immune)
TF_path="D:/project/07.STAD_subtype/analysis0507_v7/04.OS/TF"
TF=read.table(file.path(TF_path,"single_cox_sig.txt"),sep="\t",header=T,stringsAsFactors = F,quote="")
setwd("D:/project/07.STAD_subtype/analysis0507_v7/07.scRNA-seq/harmony_0.5_30/markers_by_self")
# pdf("Featureplot_TF.pdf",width=6.5,height = 10)
# FeaturePlot(data.combined.new,features = TF$Characteristics,pt.size = 0,ncol=2,reduction = "tsne")
# dev.off()
p=VlnPlot(data.combined.new, features = TF$Characteristics,pt.size = 0,stacked = T)
pdf("vlnplot_TF.pdf",width=5,height=8)
p
dev.off()
p2=DotPlot(data.combined.new,features = TF$Characteristics, cols = c("#fddbc7", "red"), dot.scale = 5) + RotatedAxis()+coord_flip()
pdf("DotPlot_TF.pdf",width=6,height=8)
p2
dev.off()
cluster_markers=read.csv("Cluster_markers.csv",header=T,stringsAsFactors = F,quote="")
cluster_markers=cluster_markers[cluster_markers$cluster=="Fibroblasts cells",]
gene_inter=intersect(TF$Characteristics,cluster_markers$gene)
#################
##----单细胞marker ssGSEA 反卷积回去----------
rm(list=ls())
data_path="D:/project/07.STAD_subtype/analysis0507_v7/00.data/mRNA"
expression = read.table(file.path(data_path,"TCGA-STAD-mRNA-tpm_log_case.txt"),header=TRUE,row.names=1, quote="",stringsAsFactors=F,check.names=F)
setwd("D:/project/07.STAD_subtype/analysis0507_v7/07.scRNA-seq/harmony_0.5_30/markers_by_self")
##读取背景基因集合
gene_set=read.csv("Cluster_markers.csv",header=T,stringsAsFactors = F)[,c(7,6)]
list=split(as.matrix(gene_set)[,1], gene_set[,2])##存储的是每个免疫细胞对应的基因,构建背景基因集合
library(GSVA)
gsva_matrix=gsva(as.matrix(expression),list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)##进行gsva分析
dir.create("ssGSEA")
write.table(gsva_matrix,"ssGSEA/gsva_matrix.txt",sep="\t",quote=FALSE,row.names = TRUE)##输出结果