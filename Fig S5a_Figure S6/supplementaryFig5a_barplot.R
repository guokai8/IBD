library(phyloseq)
library(ggplot2)
library(ggrepel)
library(VennDetail)
library(dplyr)

load(file = "makePhyloseq_percentage.rdata")
######Abundance- barplot- phylum#####
P.df<-psmelt(P.perc.ps)
ggplot(data = P.df,aes(x=ID,y=Abundance,fill=Phylum))+
  geom_bar(stat = "identity")+theme(axis.text.x = element_text(angle = 90))
mycor<-setcolor(100)
P.df%>%mutate(Group=paste(P.df$Sample,P.df$Mice,P.df$Treatment,sep = "_"))%>%
  ggplot(aes(x=Group,y=Abundance,fill=Phylum))+
  geom_bar(position = "fill",stat="identity")+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = mycor[7:50])
