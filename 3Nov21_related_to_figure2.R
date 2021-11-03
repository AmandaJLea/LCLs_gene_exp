#########
# example gene
#########

library(data.table)
library(limma)
library(qvalue)
library(mashr)

setwd("~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL")

pca_res=read.delim('31Mar21_filtered_sample_info.txt')

# sample info
info2=fread('~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2020_NovaSeq/data/MASTER_LCL_log_names.txt',stringsAsFactors=F)
info2$line_id<-paste('Line',info2$Number,sep='')
info7=read.delim('~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2019_NovaSeq/data/individuals_to_use_3Feb18.txt',stringsAsFactors=F)

info=merge(info2,info7[,c('Individual.ID','LCL.Catalog.ID.in.Coriell.system','Gender','Population')],by.y='LCL.Catalog.ID.in.Coriell.system',by.x='Name')

info$pop2<-NA
info$pop2[which(info$Population %in% c('CEU','TSI','FIN','GBR','IBS'))]<-'EUR'
info$pop2[which(info$Population %in% c('YRI','LWK','GWD','MSL','ESN','ASW','ACB'))]<-'AFR'

info=merge(info,pca_res,all.y=T,by.x='line_id',by.y='line')
info=info[order(info$id),]
info$pop2[which(is.na(info$pop2))]<-'AFR'

voom_RNA=fread('31Mar21_all_runs_voom_resid.txt')

genes=read.delim('31Mar21_all_runs_voom_norm_geneIDs.txt')
genes$ID<-1:dim(genes)[1]
genepos <- read.delim("~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/pilot_experiment/pilot_data_8Jan18/hg38_PCgene_start_end.txt", header = F)
genes2<-merge(genes,genepos,by.x='V1',by.y='V4')

subset(genes2,V5=='IL10')
tmp<-as.data.frame(t(voom_RNA[4620,]))
tmp$Population<-info$pop2
tmp$Condition<-info$treatment

library(ggplot2)
ggplot(tmp, aes(x=factor(Population), y=V1,fill=Condition)) + geom_boxplot(outlier.shape = NA) + theme_bw(13)+ scale_fill_brewer(palette="Paired")+ylim(-3,3)+theme(legend.position = "top")

tmp2<-as.data.frame(aggregate(tmp$V1~tmp$Condition+tmp$Population,FUN=mean))
names(tmp2)<-c('Condition2','Population','Expression')
tmp2$Condition<-as.character(tmp2$Condition2)
tmp2$Condition[which(tmp2$Condition=='FSL')]<-'FSL-1'
tmp2$Condition[which(tmp2$Condition=='IGF')]<-'IGF-1'

ggplot(tmp2, aes(x=factor(Population), y=Expression)) + geom_violin() +geom_boxplot(width=0.1,outlier.shape = NA)+ theme_bw(13)+ scale_color_brewer(palette="Paired")+theme(legend.position = "top")+geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.3), aes(color=factor(Condition)))

################
# ancestry effects
################

setwd("~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL")
library(mashr)
library(gplots)
library(corrplot)

lfsr=read.delim('31Mar21_LFSR_popDE_mashr.txt')
pm=read.delim('31Mar21_pm_popDE_mashr.txt')

names(pm)[5]<-'FSL-1'
names(pm)[8]<-'IGF-1'
corrplot(cor(pm),method='ellipse',type='upper')  

tmp<-apply(lfsr,2,function(x) length(which(x<0.05)))
# assess sharing
lfsr$sig<-apply(lfsr,1,function(x) length(which(x<0.1)))
lfsr$id<-1:dim(lfsr)[1]
lfsr_sig1<-subset(lfsr,sig>0 )

lfsr_sig1$within2_v2<-NA
for (i in 1:dim(lfsr_sig1)[1]){
	z<-min(lfsr_sig1[i,1:11])[1]
	focal<-pm[lfsr_sig1$id[i],which(lfsr_sig1[i,1:11] == z)][1]
	not_focal<-pm[lfsr_sig1$id[i],which(lfsr_sig1[i,1:11] != z)]
	not_focal2<-log2(not_focal/as.numeric(focal))
	lfsr_sig1$within2_v2[i]<- length(which(not_focal2 > -1 & not_focal2<1))
}

x<-as.data.frame(table(lfsr_sig1$within2_v2+1))
ggplot(x, aes(x=Var1,y=Freq )) +geom_bar(stat = "identity",fill="steelblue" )+theme_bw(13)
  
out<-as.data.frame(matrix(ncol=11,nrow=11))
rownames(out)<-colnames(out)<-names(lfsr_sig1)[1:11]

for (i in c(names(lfsr_sig1)[1:11])){
	for (k in c(names(lfsr_sig1)[1:11])){
	lfsr_tmp<-lfsr[,c(k,i)]
	pm_tmp<-pm[,c(k,i)]
	lfsr_tmp$sig<-apply(lfsr_tmp,1,function(x) length(which(x<0.1)))
	lfsr_tmp$id<-1:dim(lfsr_tmp)[1]
	lfsr_tmp2<-subset(lfsr_tmp,sig>0)
	pm_tmp$id<-1:dim(pm_tmp)[1]
	both<-merge(lfsr_tmp2,pm_tmp,by='id')
	both$log2<-log2(both[,5]/both[,6])
	out[i,k]<-dim(subset(both,log2> -1 & log2<1))[1] }}

rownames(out)[5]<-colnames(out)[5]<-'FSL-1'
rownames(out)[8]<-colnames(out)[8]<-'IGF-1'

heatmap.2(as.matrix(out),dendrogram = "none",trace="none", notecol = "Black",notecex = 2,key=T,col=colorRampPalette(c("Steelblue","White","Red")),margins=c(12,9),density.info='none')

################
# GSEA
################

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

output_coef=read.delim('31Mar21_popDE_coef.txt')

genes=read.delim('31Mar21_all_runs_voom_norm_geneIDs.txt')
protein_coding=read.delim('~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2019_NovaSeq/data/protein_coding_hg38.txt')
genes$id<-1:dim(genes)[1]
genes2=merge(genes,protein_coding,by.x='V1',by.y='Gene.stable.ID')
genes2=genes2[order(genes2$id),]

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
keytypes(org.Hs.eg.db)

output_coef$geneID<-genes2$V1
dup<-subset(as.data.frame(table(output_coef$geneID)),Freq>1)
output_coef2<-output_coef[-which(output_coef$geneID %in% dup$Var1),]
output_coef2$mean1<-apply(output_coef2[,c('BAFF','ETOH','GARD','IFNG','IGF','DEX','TUNIC','H20')],1,mean)
output_coef2$mean2<-apply(output_coef2[,c(1:12)],1,mean)

# using abs value
original_gene_list <- (output_coef2$mean1)
names(original_gene_list) <- output_coef2$geneID
gene_list<-na.omit(original_gene_list)
gene_list = sort( (gene_list), decreasing = TRUE)

gse_res <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             nPerm = 1000, 
             minGSSize = 10, maxGSSize = 500,
             pvalueCutoff = 1, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")

print(emapplot(gse_res, showCategory = 15,repel = TRUE,cex_label_category=2) ) 

gse_res2<-as.data.frame(gse_res)
write.table(gse_res2[,1:8],paste('10Sep21_GSEA_popDE_genes_mean1_NOabs.txt',sep=''),sep='\t',row.names=F)
#print(emapplot(gse_res, showCategory = 15) ) 

################
# Fst
################

cd /Genomics/ayroleslab2/alea/LCL_stim/final_dataset/

library(data.table)
library(qvalue)

geneIDs=read.delim('31Mar21_all_runs_voom_norm_geneIDs.txt')
geneIDs$id<-1:dim(geneIDs)[1]
geneIDs2=read.delim('/Genomics/ayroleslab2/alea/ref_genomes/hg38/hg38_PCgene_start_end.txt',header=F)

geneIDs_wsym<-merge(geneIDs2,geneIDs,all.y=T,by.x='V4',by.y='V1')
geneIDs_wsym<-geneIDs_wsym[order(geneIDs_wsym$id),]

data=fread('/Genomics/ayroleslab2/alea/ref_genomes/hg38/hg38_PC_genes_10kb_TSS_Fst.bed')
fst_by_gene<-as.data.frame(aggregate(data$V4~data$V9,FUN=mean))
names(fst_by_gene)<-c('gene','med_fst')
tmp<-quantile(fst_by_gene$med_fst,seq(0,1,0.01))
fst_by_gene$outlier<-0
fst_by_gene$outlier[which(fst_by_gene$med_fst>tmp[100])]<-1

# output from mashR
lfsr=read.delim('31Mar21_LFSR_popDE_mashr.txt')
pm=read.delim('31Mar21_pm_popDE_mashr.txt')

lfsr$sig<-apply(lfsr,1,function(x) length(which(x<0.1)))
lfsr$id<-1:dim(lfsr)[1]
lfsr_sig1<-lfsr

lfsr_sig1$geneID<-geneIDs_wsym$V5[lfsr_sig1$id]
lfsr_sig2<-merge(lfsr_sig1[,c('id','sig')],geneIDs_wsym,by='id',all.y=T)
lfsr_sig2<-merge(lfsr_sig2,fst_by_gene,by.x='V5',by.y='gene')
lfsr_sig2$sig2<-0
lfsr_sig2$sig2[which(lfsr_sig2$sig>1)]<-1

write.table(lfsr_sig2,'10Sep21_Fst_results.txt',row.names=F,sep='\t')

####

lfsr_sig2=read.delim('10Sep21_Fst_results.txt')
tmp<-subset(lfsr_sig2,sig==1 | sig>7)

# significant in any condition
library(ggplot2)

ggplot( data=tmp,aes( y=med_fst,x=factor(sig2))) +
    geom_violin(width=0.8) +
    geom_boxplot(width=0.1,outlier.shape = NA) +theme_bw(15)+ylim(0.025,0.4)
    
################
# Qst
################

treat1<-c('BAFF','ETOH','GARD','IFNG','IGF','DEX','TUNIC','H20')
treat2<-c('BAFF','ETOH','GARD','IFNG','IGF-1','DEX','TUNIC','H20')
par(mfrow=c(3,3))

for (k in 1:length(treat1)){
	
	null=read.delim(paste('16Jul21_background_AFR_EUR_Pst_',treat1[k],'.txt',sep=''))
	data=read.delim(paste('16Jul21_AFR_EUR_Pst_',treat1[k],'.txt',sep=''))

plot(density(null$Pst_Values),col='grey',ylim=c(0,10),lwd=3,main=treat2[k],xlab='Pst',ylab='Density',cex.lab=1.2,cex.axis=1.2)
lines(density(data$Pst_Values),col='steelblue',lwd=3)
arrows(0.1045843, 2.5, 0.1045843, 0,lwd=1,length=0.1) # 0.1045281 - 0.1046404
#legend("topleft",c("",""),col=c('steelblue','grey'),lwd=c(3,3),bty='n')
}

treat1<-c('BAFF','ETOH','GARD','IFNG','IGF','DEX','TUNIC','H20')
i=8	
	null=read.delim(paste('16Jul21_background_AFR_EUR_Pst_',treat1[k],'.txt',sep=''))
	data=read.delim(paste('16Jul21_AFR_EUR_Pst_',treat1[k],'.txt',sep=''))

par(mfrow=c(1,1))
plot(density(null$Pst_Values),col='grey',ylim=c(0,10),lwd=3,main=treat2[k],xlab='Pst',ylab='Density')
lines(density(data$Pst_Values),col='steelblue',lwd=3)
arrows(0.1045843, 2.5, 0.1045843, 0,lwd=1,length=0.1) # 0.1045281 - 0.1046404
legend("topleft",c("",""),col=c('steelblue','grey'),lwd=c(3,3),bty='n')

