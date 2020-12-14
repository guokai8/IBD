library(richR)
gene<-read.csv("IBD_loading_new.csv",row.names=1)
hsro<-buildAnnot(species="human",keytype="SYMBOL",anntype="Reactome")
sel<-gene$Approved.symbol[1:200]
res<-enrich(sel,hsro)
write.csv(res,file="reactome_enrichment.csv")
ggnetwork(res,usePadj=F)
dev.print(pdf,file="reactome_network_new.pdf")




