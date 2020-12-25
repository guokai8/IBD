library(phyloseq)
library(ggplot2)
library(ggrepel)
library(VennDetail)
library(dplyr)
library(eoffice)
library(tidyverse)
library(rstatix)
library(ggpubr)

load(file = "makePhyloseq_counts.rdata")

######alpha diversity####
Grich<-estimate_richness(G.cou.ps)
Grich$Sample_ID<-rownames(Grich)
Grich4<-left_join(Grich,as.data.frame(as.matrix(sample_data(G.cou.ps))),by=c("Sample_ID"="ID"))%>%
  select(6,7,10:12)%>%gather(index,val,-3,-4,-5)
Grich4$Mice<-factor(as.character(Grich4$Mice),levels = c("WT","Cgas"))
Grich4$Treatment<-factor(as.character(Grich4$Treatment),levels = c("Water","DSS"))
Grich4$index<-factor(as.character(Grich4$index),levels = c("Shannon","Simpson"))

Grich4%>%ggboxplot(x = "Mice",y = "val",color="Treatment",fill="Treatment",alpha=0.7)+
  scale_fill_manual(values=setcolor(2))+scale_color_manual(values=setcolor(2))+
  facet_wrap(~index,scales="free",nrow=1)+ylab(label = "Alpha Diversity")+
  theme_bw()+xlab("")+scale_x_discrete(labels=c("WT"="WT","Cgas"=bquote(Cgas^"-/-")))+
  theme(strip.background=element_blank(),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        strip.text=element_text(color="black",size = 11))+
  stat_compare_means(aes(group=Treatment,label=paste0("p=", ..p.format..)),
                     method="t.test",size=3,label.y.npc=c("bottom","top"))


Grich4%>%spread(index,val)%>%filter(Mice=="WT")%>%.[,4:5]%>%
  apply(2,function(x)t.test(x[c(1:6)],x[c(7:11)])$p.value)

Grich4%>%spread(index,val)%>%filter(Mice=="Cgas")%>%.[,4:5]%>%
  apply(2,function(x)t.test(x[c(1:5)],x[c(6:8)])$p.value)

Grich4%>%spread(index,val)%>%filter(Treatment=="Water")%>%.[,4:5]%>%
  apply(2,function(x)t.test(x[c(1:6)],x[c(7:11)])$p.value)

Grich4%>%spread(index,val)%>%filter(Treatment=="DSS")%>%.[,4:5]%>%
  apply(2,function(x)t.test(x[c(1:5)],x[c(6:8)])$p.value)
