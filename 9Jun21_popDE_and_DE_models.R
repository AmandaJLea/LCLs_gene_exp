library(data.table)
library(limma)
library(qvalue)
library(mashr)

setwd("~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL")

#########
# read in meta data and normalized counts
#########

pca_res=read.delim('31Mar21_filtered_sample_info.txt')

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

#########
# ancestry effect within each condition
#########

treat1<-c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA','DEX','TUNIC','H20')

output_coef<-matrix(ncol=12,nrow=dim(voom_RNA)[1])
output_se<-matrix(ncol=12,nrow=dim(voom_RNA)[1])
output_pval<-matrix(ncol=12,nrow=dim(voom_RNA)[1])

for (k in 1:length(treat1)){
	info2=subset(info,treatment==treat1[k] )
	mod3<-model.matrix(~info2$pop2)
	fit3 <-lmFit(voom_RNA[,info2$id,with=F],mod3)
	fit3 <- eBayes(fit3)
	print(length(which(qvalue(fit3$p.value[,2])$qvalues<0.1)))
output_coef[,k]<-fit3$coefficients[,2]
output_se[,k]<-fit3$stdev.unscaled[,2]*fit3$sigma
output_pval[,k]<-fit3$p.value[,2]
}

names(output_coef)<-c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA','DEX','TUNIC','H20')
names(output_se)<-c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA','DEX','TUNIC','H20')
names(output_pval)<-c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA','DEX','TUNIC','H20')
write.table(output_coef,'31Mar21_popDE_coef.txt',row.names=F,sep='\t')
write.table(output_se,'31Mar21_popDE_SE.txt',row.names=F,sep='\t')
write.table(output_pval,'31Mar21_popDE_pval.txt',row.names=F,sep='\t')

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
	print(length(which(qvalue(fit3$p.value[,2])$qvalues<0.1)))
output_coef[,k]<-fit3$coefficients[,2]
output_se[,k]<-fit3$stdev.unscaled[,2]*fit3$sigma
output_pval[,k]<-fit3$p.value[,2]
}

info$treatment2<-as.character(info$treatment)
info$treatment2[which(info$treatment=='ETOH')]<-'1_ETOH'

treat1<-c('DEX','TUNIC')

for (k in 1:length(treat1)){
info2=subset(info,treatment2==treat1[k] | treatment2=='1_ETOH')
	mod3<-model.matrix(~as.character(info2$treatment2)+info2$pop2)
	fit3 <-lmFit(voom_RNA[,info2$id,with=F],mod3)
	fit3 <- eBayes(fit3)
	print(length(which(qvalue(fit3$p.value[,2])$qvalues<0.1)))
output_coef[,k+10]<-fit3$coefficients[,2]
output_se[,k+10]<-fit3$stdev.unscaled[,2]*fit3$sigma
output_pval[,k+10]<-fit3$p.value[,2]
}

output_coef<-as.data.frame(output_coef[,-10])
output_se<-as.data.frame(output_se[,-10])
output_pval<-as.data.frame(output_pval[,-10])
names(output_coef)<-c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA','DEX','TUNIC')
names(output_se)<-c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA','DEX','TUNIC')
names(output_pval)<-c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA','DEX','TUNIC')
write.table(output_coef,'31Mar21_DE_coef.txt',row.names=F,sep='\t')
write.table(output_se,'31Mar21_DE_SE.txt',row.names=F,sep='\t')
write.table(output_pval,'31Mar21_DE_pval.txt',row.names=F,sep='\t')

#########
# treatment x ancestry effect within each condition
#########

treat1<-c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA')

output_coef<-matrix(ncol=12,nrow=dim(voom_RNA)[1])
output_se<-matrix(ncol=12,nrow=dim(voom_RNA)[1])

for (k in c(1:9)){
	info2=subset(info,treatment==treat1[k] | treatment=='H20')
	mod3<-model.matrix(~as.character(info2$treatment)*info2$pop2)
	fit3 <-lmFit(voom_RNA[,info2$id,with=F],mod3)
	fit3 <- eBayes(fit3)
	print(length(which(qvalue(fit3$p.value[,4])$qvalues<0.1)))
	output_coef[,k]<-fit3$coefficients[,4]
	output_se[,k]<-fit3$stdev.unscaled[,4]*fit3$sigma
}

treat1<-c('DEX','TUNIC')

for (k in 1:2){
	info2=subset(info,treatment==treat1[k] | treatment=='ETOH')
	mod3<-model.matrix(~as.character(info2$treatment)*info2$pop2)
	fit3 <-lmFit(voom_RNA[,info2$id,with=F],mod3)
	fit3 <- eBayes(fit3)
	print(length(which(qvalue(fit3$p.value[,4])$qvalues<0.1)))
	output_coef[,k+10]<-fit3$coefficients[,4]
	output_se[,k+10]<-fit3$stdev.unscaled[,4]*fit3$sigma

}

output_coef<-as.data.frame(output_coef[,-10])
output_se<-as.data.frame(output_se[,-10])
names(output_coef)<-c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA','DEX','TUNIC')
names(output_se)<-c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA','DEX','TUNIC')
write.table(output_coef,'31Mar21_popDR_coef.txt',row.names=F,sep='\t')
write.table(output_se,'31Mar21_popDR_SE.txt',row.names=F,sep='\t')

#########
# mashR - sharing of ancestry effects
#########

coef=read.delim('31Mar21_popDE_coef.txt')
se=read.delim('31Mar21_popDE_SE.txt')

# data driven + canonical covariance 
data = mash_set_data(as.matrix(coef), as.matrix(se))
U.c = cov_canonical(data)
U.pca = cov_pca(data,5)
U.ed = cov_ed(data, U.pca)
m.ed3 = mash(data, c(U.c,U.ed))

# save output
write.table(get_lfsr(m.ed3),'31Mar21_LFSR_popDE_mashr.txt',row.names=F,sep='\t')
write.table(get_pm(m.ed3),'31Mar21_pm_popDE_mashr.txt',row.names=F,sep='\t')

lfsr=read.delim('31Mar21_LFSR_popDE_mashr.txt')
pm=read.delim('31Mar21_pm_popDE_mashr.txt')

# assess sharing
lfsr$sig<-apply(lfsr,1,function(x) length(which(x<0.1)))
lfsr$id<-1:dim(lfsr)[1]
lfsr_sig1<-subset(lfsr,sig>0 )

lfsr_sig1$within2_v2<-NA
for (i in 1:dim(lfsr_sig1)[1]){
	z<-min(lfsr_sig1[i,1:11])[1]
	focal<-pm[lfsr_sig1$id[i],which(lfsr_sig1[i,1:12] == z)][1]
	not_focal<-pm[lfsr_sig1$id[i],which(lfsr_sig1[i,1:12] != z)]
	not_focal2<-log2(not_focal/as.numeric(focal))
	lfsr_sig1$within2_v2[i]<- length(which(not_focal2 > -1 & not_focal2<1))
}

barplot(table(lfsr_sig1$within2_v2+1))

#########
# mashR - ancestry x treatment effects
#########

treat1<-c('BAFF','ETOH','GARD','IFNG','IGF')

for (i in 1:length(treat1)){
	tmp<-lfsr[,which(names(lfsr)=='H20' | names(lfsr)=='id' | names(lfsr)==treat1[i])]
	tmp2<-subset(tmp,tmp[,1]<0.1 | tmp[,2]<0.1)
	tmp3<-pm[tmp2$id,which(names(lfsr)=='H20' |  names(lfsr)==treat1[i])]
	tmp3$diff<-log2(tmp3[,1]/tmp3[,2])
	print(dim(subset(tmp3,diff>1 | diff< -1)))
}

treat1<-c('DEX','TUNIC')

for (i in 1:length(treat1)){
	tmp<-lfsr[,which(names(lfsr)=='ETOH' | names(lfsr)=='id' | names(lfsr)==treat1[i])]
	tmp2<-subset(tmp,tmp[,1]<0.1 | tmp[,2]<0.1)
	tmp3<-pm[tmp2$id,which(names(lfsr)=='ETOH' |  names(lfsr)==treat1[i])]
	tmp3$diff<-log2(tmp3[,1]/tmp3[,2])
	print(dim(subset(tmp3,diff>1 | diff< -1)))
}

#########
# mashR - sharing of treatment effects
#########

coef=read.delim('31Mar21_DE_coef.txt')
se=read.delim('31Mar21_DE_SE.txt')

# data driven + canonical covariance 
# account for correlations (e.g., because conditions use the name control comparison)
data = mash_set_data(as.matrix(coef), as.matrix(se))
V = estimate_null_correlation_simple(data)
data.V = mash_update_data(data, V=V)
U.c = cov_canonical(data.V)
U.pca = cov_pca(data.V,5)
U.ed = cov_ed(data.V, U.pca)
m2.ed3 = mash(data, c(U.c,U.ed))

# save results
write.table(get_lfsr(m2.ed3),'31Mar21_LFSR_DE_mashr.txt',row.names=F,sep='\t')
write.table(get_pm(m2.ed3),'31Mar21_pm_DE_mashr.txt',row.names=F,sep='\t')

lfsr=read.delim('31Mar21_LFSR_DE_mashr.txt')
pm=read.delim('31Mar21_pm_DE_mashr.txt')

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

barplot(table(lfsr_sig1$within2_v2+1))

