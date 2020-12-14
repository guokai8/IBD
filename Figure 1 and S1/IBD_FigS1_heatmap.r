###load the package
library(DESeq2)
### load the data
## host_tx_counts.tsv was downloaded from https://ibdmdb.org/tunnel/public/HMP2/HTX/1730/products
tx<-read.table("host_tx_counts.tsv")
meta<-read.csv('hmp2_metadata-sum.csv')
### only use colon rectum and ileum samples
id<-intersect(meta$External.ID,colnames(tx))
id<-id[!id%in%c('CSMDRVXI','MSM719ME')]
dd<-tx[,id]
rownames(meta)<-meta$External.ID
condition<-meta[id,"Group"]
### build DESeq object
dds<-DESeqDataSetFromMatrix(dd,DataFrame(condition),~condition)
dds<-DESeq(dds)
vst<-varianceStabilizingTransformation(dds)
###make PCA plot
plotPCA(vst)+scale_color_manual(values=c( "#A6CEE3", "#1F78B4", "#B2DF8A")+theme_classic(base_size=15)
dev.print(pdf,file="FigureS1A_PCA.pdf")
df<-assay(vst)
group<-as.data.frame(colData(vst))
colon<-df[,rownames(subset(group,condition=="Colon"))]
rectum<-df[,rownames(subset(group,condition=="Rectum"))]
ileum<-df[,rownames(subset(group,condition=="Ileum"))]
write.table(colon,file="vst_colon.txt",sep="\t")
write.table(rectum,file="vst_rectum.txt",sep="\t")
write.table(ileum,file="vst_ileum.txt",sep="\t")
gx<-read.csv('ibd_gene.csv',header=F)
###heatmap rectum
ggr<-meta[colnames(rectum),]
ggr<-ggr[order(ggr$diagnosis),]
ann<-ggr[,"diagnosis",drop=F]
pheatmap(rectum[gx$V1,rownames(ggr)],cluster_rows=F,cluster_cols=F,scale="row",color=colorRampPalette(c('cyan4','white','red'))(256),annotation_col=ann,border_color=F,fontsize_col=1,fontsize_row=8)
dev.print(pdf,file="FigS1C_rectum.pdf")
###heatmap colon
ggr<-meta[colnames(colon),]
ggr<-ggr[order(ggr$diagnosis),]
ann<-ggr[,"diagnosis",drop=F]
pheatmap(colon[gx$V1,rownames(ggr)],cluster_rows=F,cluster_cols=F,scale="row",color=colorRampPalette(c('cyan4','white','red'))(256),annotation_col=ann,border_color=F,fontsize_col=1,fontsize_row=8)
dev.print(pdf,file="FigS1C_colon.pdf")
###heatmap ileum
ggr<-meta[colnames(ileum),]
ggr<-ggr[order(ggr$diagnosis),]
ann<-ggr[,"diagnosis",drop=F]
pheatmap(ileum[gx$V1,rownames(ggr)],cluster_rows=F,cluster_cols=F,scale="row",color=colorRampPalette(c('cyan4','white','red'))(256),annotation_col=ann,border_color=F,fontsize_col=1,fontsize_row=8)
dev.print(pdf,file="FigS1C_ileum.pdf")
