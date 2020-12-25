library(phyloseq)
library(ggplot2)

load(file ="makePhyloseq_percentage.rdata")

G.ord <- ordinate(G.perc.ps, "PCoA", "bray")
plot_ordination(G.perc.ps, G.ord,color="Treatment",shape = "Mice")+
  geom_point(size=5,alpha=.7)+geom_text_repel(aes(label=sample_names(G.perc.ps)),color="black",size=2.5)
G.pcoa<-plot_ordination(G.perc.ps, G.ord,color="Treatment",shape = "Mice")

G.pcoa$data%>%ggplot(aes(x=Axis.1,y=Axis.2,color=Treatment,shape=Mice))+
  scale_shape_manual(values = c(16,17))+geom_point(size=4,alpha=0.8,)+
  stat_ellipse()+scale_color_manual(values = c("#377eb8","#e41a1c"))+
  theme_light(base_size = 15)+xlab("PC1 (34.4%)")+ylab("PC2 (20.3%)")

library(vegan)
set.seed(2)
adonis(phyloseq::distance(G.perc.ps, method = "bray") ~ Mice+Treatment, data = data.frame(sample_data(G.perc.ps)))

