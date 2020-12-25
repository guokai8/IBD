library(phyloseq)
library(ggplot2)
library(ggrepel)
library(VennDetail)
library(dplyr)
library(tidyverse)
library(rstatix)
library(ggpubr)

load(file = "makePhyloseq_counts.rdata")
#####DEseq#######
library(DESeq2)
G_WT_d2<-phyloseq_to_deseq2(subset_samples(G.cou.ps, Mice=="WT"),~Treatment)
gm_mean=function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]),na.rm=na.rm)/length(x))
}             
G_WT_d2=estimateSizeFactors(G_WT_d2,geoMeans=apply(counts(G_WT_d2), 1, gm_mean))%>%DESeq(.,fitType="local")  #normalize#
G_WT_d2=results(G_WT_d2)%>%.[order(.$pvalue, na.last=NA), ]
G_WT_rr=cbind(as.data.frame(G_WT_d2),as.matrix(tax_table(G.cou.ps)[rownames(G_WT_d2),]))
fun_rr<-function(rr){
  sig_otu<-rr%>%rownames_to_column(var="sig_otu")%>%filter(pvalue<0.05,abs(log2FoldChange)>1)%>%pull(sig_otu)
  rr[sig_otu,"sig_otu"]<-sig_otu
  rr$label<-paste(rr$Genus," (",rr$sig_otu,")",sep="")
  table(rr$sig_otu=="",exclude = NULL)
  rr$color<-0
  rr[rr%>%filter(pvalue<0.05,log2FoldChange<0)%>%pull(sig_otu),"color"]<-2
  rr[rr%>%filter(pvalue<0.05,log2FoldChange>0)%>%pull(sig_otu),"color"]<-1
  return(rr)
}
G_WT_rr<-fun_rr(G_WT_rr) 
G_WT_rr$label[18:nrow(G_WT_rr)]<-""
#supplementary fig 5b
G_WT_rr%>%filter(color!=0)%>%ggplot(aes(x = log2FoldChange,y = Genus,color=factor(color)))+
  geom_vline(xintercept=0.0,color="gray",size=0.5)+geom_point(size=3,alpha=0.8)+
  theme(axis.text.x=element_text(angle=0,hjust = 0,vjust=0.5))+theme_bw()+
  scale_color_manual(values = c("#e41a1c","#377eb8"),labels=c("Up","Down"))+
  guides(color=guide_legend(title=NULL))+ggtitle(label = "WT: DSS vs Water")+
  scale_y_discrete(limits=c(G_WT_rr%>%filter(color!=0)%>%.$Genus%>%sort(decreasing = T)))  ##reorder y axis
# write.table(x = G_WT_rr[,c(-9,-10)],file = "DifferentialGenus_WT.csv",quote = F,sep = "\t",col.names = T,row.names = T)
#supplementary fig 5c
G_Cgas_d2<-phyloseq_to_deseq2(subset_samples(G.cou.ps, Mice=="Cgas"),~Treatment)
G_Cgas_d2=estimateSizeFactors(G_Cgas_d2,geoMeans=apply(counts(G_Cgas_d2), 1, gm_mean))%>%DESeq(.,fitType="local")
G_Cgas_d2=results(G_Cgas_d2)%>%.[order(.$padj, na.last=NA), ]
G_Cgas_rr=cbind(as.data.frame(G_Cgas_d2),as.matrix(tax_table(G.cou.ps)[rownames(G_Cgas_d2),]))
G_Cgas_rr<-fun_rr(G_Cgas_rr) 
G_Cgas_rr$label[17:nrow(G_Cgas_rr)]<-""
G_Cgas_rr%>%filter(color!=0)%>%ggplot(aes(x = log2FoldChange,y = Genus,color=factor(color)))+
  geom_vline(xintercept=0.0,color="gray",size=0.5)+geom_point(size=3,alpha=0.8)+
  theme(axis.text.x=element_text(angle=0,hjust = 0,vjust=0.5))+theme_bw()+
  scale_color_manual(values = c("#e41a1c","#377eb8"),labels=c("Up","Down"))+
  guides(color=guide_legend(title=NULL))+ggtitle(label = "Cgas: DSS vs Water")+
  scale_y_discrete(limits=c(G_Cgas_rr%>%filter(color!=0)%>%.$Genus%>%sort(decreasing = T)))  
# write.table(x = G_Cgas_rr[,c(-9,-10)],file = "DifferentialGenus_Cgas.csv",quote = F,sep = "\t",col.names = T,row.names = T)
#figure 3d
G_Water_d2<-phyloseq_to_deseq2(subset_samples(G.cou.ps, Treatment=="Water"),~Mice)
G_Water_d2=estimateSizeFactors(G_Water_d2,geoMeans=apply(counts(G_Water_d2), 1, gm_mean))%>%DESeq(.,fitType="local")
G_Water_d2=results(G_Water_d2)%>%.[order(.$padj, na.last=NA), ]
G_Water_rr=cbind(as.data.frame(G_Water_d2),as.matrix(tax_table(G.cou.ps)[rownames(G_Water_d2),]))
G_Water_rr<-fun_rr(G_Water_rr) 
G_Water_rr$label[16:nrow(G_Water_rr)]<-""
G_Water_rr%>%filter(color!=0)%>%ggplot(aes(x = log2FoldChange,y = Genus,color=factor(color)))+
  geom_vline(xintercept=0.0,color="gray",size=0.5)+geom_point(size=3,alpha=0.8)+
  theme(axis.text.x=element_text(angle=0,hjust = 0,vjust=0.5))+theme_bw()+
  scale_color_manual(values = c("#e41a1c","#377eb8"),labels=c("Up","Down"))+
  guides(color=guide_legend(title=NULL))+ggtitle(label = "Water: Cgas vs WT")+
  scale_y_discrete(limits=c(G_Water_rr%>%filter(color!=0)%>%.$Genus%>%sort(decreasing = T)))  
# write.table(x = G_Water_rr[,c(-9,-10)],file = "DifferentialGenus_Water.csv",quote = F,sep = "\t",col.names = T,row.names = T)
#figure 3e
G_DSS_d2<-phyloseq_to_deseq2(subset_samples(G.cou.ps, Treatment=="DSS"),~Mice)
G_DSS_d2=estimateSizeFactors(G_DSS_d2,geoMeans=apply(counts(G_DSS_d2), 1, gm_mean))%>%DESeq(.,fitType="local")
G_DSS_d2=results(G_DSS_d2)%>%.[order(.$padj, na.last=NA), ]
G_DSS_rr=cbind(as.data.frame(G_DSS_d2),as.matrix(tax_table(G.cou.ps)[rownames(G_DSS_d2),]))
G_DSS_rr<-fun_rr(G_DSS_rr) 
G_DSS_rr$label[19:nrow(G_DSS_rr)]<-""
G_DSS_rr%>%filter(color!=0)%>%ggplot(aes(x = log2FoldChange,y = Genus,color=factor(color)))+
  geom_vline(xintercept=0.0,color="gray",size=0.5)+geom_point(size=3,alpha=0.8)+
  theme(axis.text.x=element_text(angle=0,hjust = 0,vjust=0.5))+theme_bw()+
  scale_color_manual(values = c("#e41a1c","#377eb8"),labels=c("Up","Down"))+
  guides(color=guide_legend(title=NULL))+ggtitle(label = "DSS: Cgas vs WT")+
  scale_y_discrete(limits=c(G_DSS_rr%>%filter(color!=0)%>%.$Genus%>%sort(decreasing = T)))  
# write.table(x = G_DSS_rr[,c(-9,-10)],file = "DifferentialGenus_DSS.csv",quote = F,sep = "\t",col.names = T,row.names = T)

