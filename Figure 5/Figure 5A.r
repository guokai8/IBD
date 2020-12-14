##load data
cgas<-read.delim('cgas.txt',sep="\t")
wt<-read.csv('wtcgas.txt',sep="\t")
##top and bottom 100
cc<-rbind(head(cgas,100),tail(cgas,100))
ww<-rbind(head(wt,100),tail(wt,100))
cc<-cc[,c('Name','Score')]
ww<-ww[,c('Name','Score')]
library(VennDetail)
ven<-venndetail(list('WT Score'=ww$Name,'Cgas-/- Score'=cc$Name))
rownames(cc)<-cc$Name
rownames(ww)<-ww$Name
##make the Venn Diagram figure
plot(ven)
dev.print(pdf,file="Figure5A_venn.pdf")
##get the Shared durg
res<-getFeature(ven,subset="Shared",rlist=list('WT Score'=ww,'Cgas-/- Score'=cc))
res<-res[,c(2,4,6)]
rownames(res)<-res[,1]
library(pheatmap)
colnames(res)[2:3]<-sub('_Score','',colnames(res)[2:3])
pheatmap(res[,2:3],cluster_cols=F,colorRampPalette(c("darkcyan","white","red"))(256),border_color = "white",display_numbers=T)
dev.print(pdf,file="Figure5A.pdf")
