library(data.table)
library(limma)
library(qvalue)
library(mashr)

setwd("~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL")

#########
# read in meta data and normalized counts
#########

# sample info
info2=fread('~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2020_NovaSeq/data/MASTER_LCL_log_names.txt')
info2$line_id<-paste('Line',info2$Number,sep='')
info7=read.delim('~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2019_NovaSeq/data/individuals_to_use_3Feb18.txt')

info=merge(info2,info7[,c('Individual.ID','LCL.Catalog.ID.in.Coriell.system','Gender','Population')],by.y='LCL.Catalog.ID.in.Coriell.system',by.x='Name')

info$pop2<-NA
info$pop2[which(info$Population %in% c('CEU','TSI','FIN','GBR','IBS'))]<-'EUR'
info$pop2[which(info$Population %in% c('YRI','LWK','GWD','MSL','ESN','ASW','ACB'))]<-'AFR'

info=merge(info,pca_res,all.y=T,by.x='line_id',by.y='line')
info=info[order(info$id),]
info$pop2[which(is.na(info$pop2))]<-'AFR'

voom_RNA=fread('31Mar21_all_runs_voom_resid.txt')
genes=read.delim('31Mar21_all_runs_voom_norm_geneIDs.txt')

#########
# ancestry effect within each condition
#########

treat1<-c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA','DEX','TUNIC','H20')

for (k in 1:length(treat1)){
	info2=subset(info,treatment==treat1[k] )
	mod3<-model.matrix(~info2$pop2)
	fit3 <-lmFit(voom_RNA[,info2$id,with=F],mod3)
	fit3 <- eBayes(fit3)
	top<-as.data.frame(topTable(fit3,number=dim(voom_RNA)[1]))
	top$gene<-as.numeric(rownames(top))
	top2<-top[order((top$gene)),]
	top2$SE<-fit3$stdev.unscaled[,2]*fit3$sigma
	top2$gene<-genes$V1
	
write.table(top2[,c(1,3,6,7)],paste('23May22_limma_ancestry_',treat1[k],'.txt',sep=''),row.names=F,sep='\t',quote=F)
}

#########
# treatment effect within each condition
#########

info$treatment2<-as.character(info$treatment)
info$treatment2[which(info$treatment=='H20')]<-'1_H20'

treat1<-c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA')

output_coef<-matrix(ncol=12,nrow=dim(voom_RNA)[1])
output_se<-matrix(ncol=12,nrow=dim(voom_RNA)[1])
output_pval<-matrix(ncol=12,nrow=dim(voom_RNA)[1])

for (k in 1:length(treat1)){
	info2=subset(info,treatment2==treat1[k] | treatment2=='1_H20')
	mod3<-model.matrix(~as.character(info2$treatment2)+info2$pop2)
	fit3 <-lmFit(voom_RNA[,info2$id,with=F],mod3)
	fit3 <- eBayes(fit3)
	top<-as.data.frame(topTable(fit3,number=dim(voom_RNA)[1]))
	top$gene<-as.numeric(rownames(top))
	top2<-top[order((top$gene)),]
	names(top2)[1]<-'Treatment_logFC'
	names(top2)[2]<-'Ancestry_logFC'
	top2$SE<-fit3$stdev.unscaled[,2]*fit3$sigma
	top2$gene<-genes$V1
write.table(top2[,c(1,2,4,6,7)],paste('23May22_limma_treatment_',treat1[k],'.txt',sep=''),row.names=F,sep='\t',quote=F)

}

info$treatment2<-as.character(info$treatment)
info$treatment2[which(info$treatment=='ETOH')]<-'1_ETOH'

treat1<-c('DEX','TUNIC')

for (k in 1:length(treat1)){
info2=subset(info,treatment2==treat1[k] | treatment2=='1_ETOH')
	mod3<-model.matrix(~as.character(info2$treatment2)+info2$pop2)
	fit3 <-lmFit(voom_RNA[,info2$id,with=F],mod3)
	fit3 <- eBayes(fit3)
	top<-as.data.frame(topTable(fit3,number=dim(voom_RNA)[1]))
	top$gene<-as.numeric(rownames(top))
	top2<-top[order((top$gene)),]
	names(top2)[1]<-'Treatment_logFC'
	names(top2)[2]<-'Ancestry_logFC'
	top2$SE<-fit3$stdev.unscaled[,2]*fit3$sigma
	top2$gene<-genes$V1
write.table(top2[,c(1,2,4,6,7)],paste('23May22_limma_treatment_',treat1[k],'.txt',sep=''),row.names=F,sep='\t',quote=F)
}

