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
S_WT_d2<-phyloseq_to_deseq2(subset_samples(S.cou.ps, Mice=="WT"),~Treatment)
gm_mean=function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]),na.rm=na.rm)/length(x))
}             
S_WT_d2=estimateSizeFactors(S_WT_d2,geoMeans=apply(counts(S_WT_d2), 1, gm_mean))%>%DESeq(.,fitType="local")  #normalize#
S_WT_d2=results(S_WT_d2)%>%.[order(.$pvalue, na.last=NA), ]
S_WT_rr=cbind(as.data.frame(S_WT_d2),as.matrix(tax_table(S.cou.ps)[rownames(S_WT_d2),]))
fun_rr<-function(rr){
  sig_otu<-rr%>%rownames_to_column(var="sig_otu")%>%filter(pvalue<0.05,abs(log2FoldChange)>1)%>%pull(sig_otu)
  rr[sig_otu,"sig_otu"]<-sig_otu
  rr$label<-paste(rr$Species," (",rr$sig_otu,")",sep="")
  table(rr$sig_otu=="",exclude = NULL)
  rr$color<-0
  rr[rr%>%filter(sig_otu!="NA",pvalue<0.05,log2FoldChange<0)%>%pull(sig_otu),"color"]<-2
  rr[rr%>%filter(sig_otu!="NA",pvalue<0.05,log2FoldChange>0)%>%pull(sig_otu),"color"]<-1
  return(rr)
}
S_WT_rr<-fun_rr(S_WT_rr) 
#supplementaryFig 6b
S_WT_rr%>%filter(padj<0.05,color!=0)%>%ggplot(aes(x = log2FoldChange,y = Species,color=factor(color)))+
  geom_vline(xintercept=0.0,color="gray",size=0.5)+geom_point(size=3,alpha=0.8)+
  theme(axis.text.x=element_text(angle=0,hjust = 0,vjust=0.5))+theme_bw()+
  scale_color_manual(values = c("#e41a1c","#377eb8"),labels=c("Up","Down"))+
  guides(color=guide_legend(title=NULL))+ggtitle(label = "WT: DSS vs Water")+
  scale_y_discrete(limits=c(S_WT_rr%>%filter(padj<0.05,color!=0)%>%.$Species%>%sort(decreasing = T)))  ##reorder y axis
# write.csv(x = S_WT_rr[,c(-9,-10)],file = "Based on Counts(Species)/DifferentialSpecies_WT.csv",quote = F)

#supplementaryFig 6c
S_Cgas_d2<-phyloseq_to_deseq2(subset_samples(S.cou.ps, Mice=="Cgas"),~Treatment)
S_Cgas_d2=estimateSizeFactors(S_Cgas_d2,geoMeans=apply(counts(S_Cgas_d2), 1, gm_mean))%>%DESeq(.,fitType="local")
S_Cgas_d2=results(S_Cgas_d2)%>%.[order(.$padj, na.last=NA), ]
S_Cgas_rr=cbind(as.data.frame(S_Cgas_d2),as.matrix(tax_table(S.cou.ps)[rownames(S_Cgas_d2),]))
S_Cgas_rr<-fun_rr(S_Cgas_rr) 
S_Cgas_rr%>%filter(padj<0.05,color!=0)%>%ggplot(aes(x = log2FoldChange,y = Species,color=factor(color)))+
  geom_vline(xintercept=0.0,color="gray",size=0.5)+geom_point(size=3,alpha=0.8)+
  theme(axis.text.x=element_text(angle=0,hjust = 0,vjust=0.5))+theme_bw()+
  scale_color_manual(values = c("#e41a1c","#377eb8"),labels=c("Up","Down"))+
  guides(color=guide_legend(title=NULL))+ggtitle(label = "Cgas: DSS vs Water")+
  scale_y_discrete(limits=c(S_Cgas_rr%>%filter(padj<0.05,color!=0)%>%.$Species%>%sort(decreasing = T)))  
# write.csv(x = S_Cgas_rr[,c(-9,-10)],file = "Based on Counts(Species)/DifferentialSpecies_Cgas.csv",quote = F)

#supplementaryFig 6a
S_Water_d2<-phyloseq_to_deseq2(subset_samples(S.cou.ps, Treatment=="Water"),~Mice)
S_Water_d2=estimateSizeFactors(S_Water_d2,geoMeans=apply(counts(S_Water_d2), 1, gm_mean))%>%DESeq(.,fitType="local")
S_Water_d2=results(S_Water_d2)%>%.[order(.$padj, na.last=NA), ]
S_Water_rr=cbind(as.data.frame(S_Water_d2),as.matrix(tax_table(S.cou.ps)[rownames(S_Water_d2),]))
S_Water_rr<-fun_rr(S_Water_rr) 
S_Water_rr%>%filter(padj<0.05,color!=0)%>%ggplot(aes(x = log2FoldChange,y = Species,color=factor(color)))+
  geom_vline(xintercept=0.0,color="gray",size=0.5)+geom_point(size=3,alpha=0.8)+
  theme(axis.text.x=element_text(angle=0,hjust = 0,vjust=0.5))+theme_bw()+
  scale_color_manual(values = c("#e41a1c","#377eb8"),labels=c("Up","Down"))+
  guides(color=guide_legend(title=NULL))+ggtitle(label = "Water: Cgas vs WT")+
  scale_y_discrete(limits=c(S_Water_rr%>%filter(padj<0.05,color!=0)%>%.$Species%>%sort(decreasing = T)))  
# write.csv(x = S_Water_rr[,c(-9,-10)],file = "Based on Counts(Species)/DifferentialSpecies_Water.csv",quote = F)

#supplementaryFig 6d
S_DSS_d2<-phyloseq_to_deseq2(subset_samples(S.cou.ps, Treatment=="DSS"),~Mice)
S_DSS_d2=estimateSizeFactors(S_DSS_d2,geoMeans=apply(counts(S_DSS_d2), 1, gm_mean))%>%DESeq(.,fitType="local")
S_DSS_d2=results(S_DSS_d2)%>%.[order(.$padj, na.last=NA), ]
S_DSS_rr=cbind(as.data.frame(S_DSS_d2),as.matrix(tax_table(S.cou.ps)[rownames(S_DSS_d2),]))
S_DSS_rr<-fun_rr(S_DSS_rr) 
S_DSS_rr%>%filter(padj<0.05,color!=0)%>%ggplot(aes(x = log2FoldChange,y = Species,color=factor(color)))+
  geom_vline(xintercept=0.0,color="gray",size=0.5)+geom_point(size=3,alpha=0.8)+
  theme(axis.text.x=element_text(angle=0,hjust = 0,vjust=0.5))+theme_bw()+
  scale_color_manual(values = c("#e41a1c","#377eb8"),labels=c("Up","Down"))+
  guides(color=guide_legend(title=NULL))+ggtitle(label = "Cgas: DSS vs Water")+
  scale_y_discrete(limits=c(S_DSS_rr%>%filter(padj<0.05,color!=0)%>%.$Species%>%sort(decreasing = T)))  
# write.csv(x = S_DSS_rr[,c(-9,-10)],file = "Based on Counts(Species)/DifferentialSpecies_DSS.csv",quote = F)

