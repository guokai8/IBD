###load data
library(tidyverse)
library(DESeq2)
library(phyloseq)
## host_tx_counts.tsv was downloaded from https://ibdmdb.org/tunnel/public/HMP2/HTX/1730/products
tx<-read.table("host_tx_counts.tsv")
ph<-phyloseq::import_biom("taxonomic_profiles.biom")
####load sample information
## hmp2_metadata.csv was downloaded from https://ibdmdb.org/tunnel/public/summary.html
meta<-read.csv("hmp2_metadata.csv")
###Only focus Rectum
sample_micro<-meta%>%filter(biopsy_location=="Rectum",data_type=="biopsy_16S"|data_type=="stool_16S",diagnosis!="CD")%>%dplyr::select(External.ID,Participant.ID,site_sub_coll,diagnosis)
sample_rna<-meta%>%filter(biopsy_location=="Rectum",data_type=="host_transcriptomics",diagnosis!="CD")%>%dplyr::select(External.ID,Participant.ID,site_sub_coll,diagnosis)
###Only use the unique ID
sample_rna<-sample_rna%>%group_by(Participant.ID)%>%dplyr::slice(1)
###Only use the unique ID
sample_micro<-sample_micro%>%group_by(Participant.ID)%>%dplyr::slice(1)
### Find the shared sample with RNA-seq and microbiome data
over<-intersect(sample_rna$Participant.ID,sample_micro$Participant.ID)
sample_rna<-sample_rna%>%filter(Participant.ID%in%over)
sample_micro<-sample_micro%>%filter(Participant.ID%in%over)
sample_rna<-as.data.frame(sample_rna)
sample_micro<-as.data.frame(sample_micro)
rownames(sample_micro)<-sample_micro$External.ID
### 
rna<-tx[,intersect(sample_rna$External.ID,colnames(tx))]
otu_table(ph)<-otu_table(ph)[,sample_micro$External.ID]
colnames(otu_table(ph))<-sample_micro[colnames(otu_table(ph)),"Participant.ID"]
rownames(sample_micro)<-sample_micro$Participant.ID
###Just keep the ID and group
sample_data(ph)<-sample_micro[,c(2,4)]
rownames(sample_rna)<-sample_rna$External.ID
colnames(rna)<-sample_rna[colnames(rna),"Participant.ID"]
rownames(sample_rna)<-sample_rna$Participant.ID
condition<-sample_rna$diagnosis[sample_rna$Participant.ID!="M2064"]
####the "M2064" contain all zero
rna<-rna[,colnames(rna)!="M2064"]
ph<-subset_samples(ph,Participant.ID!="M2064")
###find the differential expressed genes
dds<-DESeqDataSetFromMatrix(rna,DataFrame(condition),~condition)
dds<-DESeq2::DESeq(dds)
res<-results(dds,contrast = c("condition","UC","nonIBD"))
rld<-rlogTransformation(dds)
####### select DEGs and get the expression values
gene<-assay(rld)[rownames(subset(res,padj<0.01)),]
##normalize the microbiome data
phn<-transform_sample_counts(ph,function(x)x/sum(x))
micro<-as.data.frame(otu_table(phn))
micro<-micro[apply(micro, 1, function(x)sum(x>0)>5),]
###Prepare data for O2PLS
gene<-t(gene)
gene<-scale(gene)
micro<-t(micro)
micro<-scale(micro)
###load the package
library(OmicsPLS)
set.seed(123)
##determine the best paramaters
crossval_o2m_adjR2(gene,micro,1:4,0:4,0:4,nr_cores = 30,nr_folds = 10)
fit<-o2m(gene,micro,n=1,nx=2,ny=0)
###write out the data
write.csv(fit$W.,"IBD_gene_loading.csv")
write.csv(fit$C.,"IBD_microbiome_loading.csv")
###Supplemental Figure 1B
sc<-as.data.frame(cbind(fit$Tt,fit$U))
sc$Group<-condition
p<-ggplot(sc,aes(V1,V2,color=Group))+geom_point()+xlab("Joint transcript score")+ylab('Joint microbiome score')+theme_light(base_size=15)
ggsave(p,file="Supplemental Figure 1B.pdf")
####
