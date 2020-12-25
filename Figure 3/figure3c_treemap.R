######Abundance- treeplot- phylum####
library(treemapify)
library(phyloseq)
library(ggplot2)
library(ggrepel)
library(VennDetail)
library(dplyr)

load(file ="makePhyloseq_percentage.rdata")

P.df<-psmelt(P.perc.ps)
P.df$phylum<-factor(P.df$Phylum,levels = c(P.df%>%arrange(desc(Abundance))%>%.$Phylum%>%unique())) #####order Phylum by Abundance###
P.df[P.df$Abundance>0.1,]%>%
  group_by(Mice,Treatment,phylum)%>%summarise(Sum=sum(Abundance))%>%
  mutate(Prec=paste(phylum,"\n",round(100*Sum/sum(Sum),digits = 2),"%"))%>%
  ggplot(aes(area=Sum,fill=phylum,label=Prec))+geom_treemap()+
  facet_grid(Mice~Treatment)+labs(fill = "Phylum")+
  scale_fill_manual(values = c("#8dd3c7","#b3de69","#fccde5","#ADE8F4","#fdb462","#74C69D",
                               "#F7D9C4","#80b1d3","#FF928B","#bc80bd","#d9d9d9","#8dd3c7"))+
  geom_treemap_text(colour = "black", place = "centre",grow = TRUE)+
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size=12,face = "bold"),legend.title = element_text(size=12))
# ggsave("PhylumAbundance-treemap.pdf",width = 6,height = 4)
