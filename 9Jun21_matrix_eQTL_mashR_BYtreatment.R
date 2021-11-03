library(data.table)
library(qvalue)

treatments<-c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA','DEX','TUNIC','H20')

for (k in 1:length(treatments)){
j=1

data=fread(paste("31Mar21_LCLs_cis_eQTL_chr",j,".",treatments[k],"0.matEQTL",sep=''))
for (j in 2:22){
data_tmp=fread(paste("31Mar21_LCLs_cis_eQTL_chr",j,".",treatments[k],"0.matEQTL",sep=''))
data=rbind(data,data_tmp) }

data$treatment<-treatments[k]
names(data)[5]<-'p_value'
data$FDR2<-qvalue(data$p_value)$qvalues
data$sig<-0
data$sig[which(data$FDR2<0.2)]<-1

write.table(data,paste("31Mar21_LCLs_cis_eQTL_ALLchr",".",treatments[k],"0.matEQTL",sep=''))
print(k)}

### get all SNPs tested

k-1
j=1

data=fread(paste("31Mar21_LCLs_cis_eQTL_chr",j,".",treatments[k],"0.matEQTL",sep=''))
for (j in 2:22){
data_tmp=fread(paste("31Mar21_LCLs_cis_eQTL_chr",j,".",treatments[k],"0.matEQTL",sep=''))
data=rbind(data,data_tmp) }

snps=unique(data[,c('SNP','gene'),with=F])
write.table(snps,'31Mar21_LCLs_all_tested_SNPs.txt',row.names=F,sep='\t')
###

k=12
data=fread(paste("31Mar21_LCLs_cis_eQTL_ALLchr",".",treatments[k],"0.matEQTL",sep=''))
data$treatment<-treatments[k]
rand<-sample(1:dim(data)[1],50000)

for (k in 1:11){
tmp_data=fread(paste("31Mar21_LCLs_cis_eQTL_ALLchr",".",treatments[k],"0.matEQTL",sep=''))
tmp_data$treatment<-treatments[k]
data=rbind(data,tmp_data)
}

data$gene_SNP<-paste(data$gene,data$SNP,sep='_')
tmp<-subset(data,sig==1)
strong_sites<-as.data.frame(table(tmp$gene_SNP))
strong_sites2<-subset(data,gene_SNP %in% strong_sites$Var1)

library(reshape2)
strong_sites2$se<-strong_sites2[,4]/strong_sites2[,5]
data_wide <- dcast(strong_sites2, SNP +gene ~ treatment, value.var="beta")
data_wide2 <- dcast(strong_sites2, SNP +gene ~ treatment, value.var="se")
data_wide3 <- dcast(strong_sites2, SNP +gene ~ treatment, value.var="FDR2")

write.table(data_wide,"31Mar21_LCLs_strong_eQTL_ALLchr_beta.txt",row.names=F,sep='\t')
write.table(data_wide2,"31Mar21_LCLs_strong_eQTL_ALLchr_SE.txt",row.names=F,sep='\t')
write.table(data_wide3,"31Mar21_LCLs_strong_eQTL_ALLchr_FDR.txt",row.names=F,sep='\t')

tmp<-subset(data,treatment=='DEX')$gene_SNP[rand]
rand_sites<-subset(data,gene_SNP %in% tmp)

library(reshape2)
rand_sites$se<-rand_sites[,4]/rand_sites[,5]
data_wide_rand <- dcast(rand_sites, SNP +gene ~ treatment, value.var="beta")
data_wide2_rand <- dcast(rand_sites, SNP +gene ~ treatment, value.var="se")
data_wide3_rand <- dcast(rand_sites, SNP +gene ~ treatment, value.var="FDR2")

write.table(data_wide_rand,"31Mar21_LCLs_rand_eQTL_ALLchr_beta.txt",row.names=F,sep='\t')
write.table(data_wide2_rand,"31Mar21_LCLs_rand_eQTL_ALLchr_SE.txt",row.names=F,sep='\t')
write.table(data_wide3_rand,"31Mar21_LCLs_rand_eQTL_ALLchr_FDR.txt",row.names=F,sep='\t')

# how many eQTL genes from matrix eQTL?
data=read.delim("31Mar21_LCLs_strong_eQTL_ALLchr_beta.txt")
for (i in 3:14){
tmp<-subset(data,data[,i]<0.1)
# print(length(unique(tmp$SNP)))
print(length(unique(tmp$gene)))
}

##########
# mashR
##########

# https://stephenslab.github.io/mashr/articles/eQTL_outline.html
library(gplots)
library(ggplot2)
library(mashr)
library(data.table)

data_wide=fread("31Mar21_LCLs_strong_eQTL_ALLchr_beta.txt")
data_wide2=fread("31Mar21_LCLs_strong_eQTL_ALLchr_SE.txt")
data_wide3=fread("31Mar21_LCLs_strong_eQTL_ALLchr_FDR.txt")

data_wide_rand=fread("31Mar21_LCLs_rand_eQTL_ALLchr_beta.txt")
data_wide2_rand=fread("31Mar21_LCLs_rand_eQTL_ALLchr_SE.txt")
data_wide3_rand=fread("31Mar21_LCLs_rand_eQTL_ALLchr_FDR.txt")

U.pca = cov_pca(data.strong,5)
U.ed = cov_ed(data.strong, U.pca)
U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)

write.table(get_lfsr(m2),'31Mar21_LFSR_eQTL_mashr_nocor.txt',row.names=F,sep='\t')
write.table(get_pm(m2),'31Mar21_pm_eQTL_mashr_nocor.txt',row.names=F,sep='\t')

lfsr=read.delim('31Mar21_LFSR_eQTL_mashr_nocor.txt')
pm=read.delim('31Mar21_pm_eQTL_mashr_nocor.txt')

lfsr$sig<-apply(lfsr,1,function(x) length(which(x<0.1)))
lfsr$id<-1:dim(lfsr)[1]
lfsr$gene<-data$gene
lfsr$SNP<-data$SNP

# how many eQTL genes?
for (i in 1:12){
tmp<-subset(lfsr,lfsr[,i]<0.1)
# print(length(unique(tmp$SNP)))
print(length(unique(tmp$gene)))
}

# sig in at least one condition from mashR
lfsr$mean<-apply(lfsr[,1:12],1,function(x) length(unique(x)))
lfsr_sig1<-subset(lfsr,sig>0 & mean>1)
lfsr_sig1$within2_v2<-NA
lfsr_sig1$focal_pm<-NA
lfsr_sig1$other_pm<-NA

for (i in 1:dim(lfsr_sig1)[1]){
	z<-min(lfsr_sig1[i,1:12])[1]
	focal<-pm[lfsr_sig1$id[i],which(lfsr_sig1[i,1:12] == z)][1]
	not_focal<-pm[lfsr_sig1$id[i],which(lfsr_sig1[i,1:12] != z)]
	not_focal2<-log2(not_focal/as.numeric(focal))
	lfsr_sig1$within2_v2[i]<-length(which(not_focal2 > -1 & not_focal2<1))
	lfsr_sig1$focal_pm[i]<-mean(focal)
	lfsr_sig1$other_pm[i]<-mean(t(not_focal))
}

table(lfsr_sig1$within2_v2+1)
tmp2<-unique(lfsr_sig1[,c('gene','within2_v2')])
table(tmp2$within2_v2)

write.table(as.data.frame(lfsr_sig1),'10Jun21_eQTL_SNPs_sharing.txt',row.names=F,sep='\t')

# not shared with other conditions
lfsr_sig2<-subset(lfsr_sig1,sig==1 & within2_v2==0)
lfsr_sig2$treatment<-apply(lfsr_sig2[,1:12],1,function(x) names(lfsr_sig2)[which(x==min(x))] ) 
as.matrix(table(lfsr_sig2$treatment))

tmp<-unique(lfsr_sig2[,c('treatment','gene')])
as.matrix(table(tmp$treatment))

###########
# response eQTL
###########

data_wide3=read.delim("31Mar21_LCLs_strong_eQTL_ALLchr_FDR.txt")
lfsr=read.delim('31Mar21_LFSR_eQTL_mashr_nocor.txt')
pm=read.delim('31Mar21_pm_eQTL_mashr_nocor.txt')

lfsr$gene<-data_wide3$gene
lfsr$SNP<-data_wide3$SNP
pm$gene<-data_wide3$gene
pm$SNP<-data_wide3$SNP
lfsr$sig<-apply(lfsr,1,function(x) length(which(x<0.1)))
lfsr$id<-1:dim(lfsr)[1]
pm$id<-1:dim(pm)[1]

out_diff<-c()
out_id<-c()
out_treat<-c()
out_pm_treat<-c()
out_pm_control<-c()

# treatment vs control
treatments<-c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA')

for (i in 1:length(treatments)){
	tmp<-lfsr[,which(names(lfsr)==treatments[i] | names(lfsr)=='H20' | names(lfsr)=='id')]
	tmp2<-subset(tmp,tmp[,1]<0.1 | tmp[,2]<0.1)
	tmp3<-pm[tmp2$id,which(names(pm)==treatments[i] | names(pm)=='H20' | names(pm)=='id' )]
	tmp3$diff<-log2(tmp3[,1]/tmp3[,2])
	tmp4<-subset(tmp3,diff < -1 | diff >1)
	out_diff<-c(out_diff,tmp4$diff)
	out_id<-c(out_id,tmp4$id)
	out_treat<-c(out_treat,as.character(rep(treatments[i],dim(tmp4)[1])))
	out_pm_treat<-c(out_pm_treat,tmp4[,which(names(tmp4)==treatments[i])])
	out_pm_control<-c(out_pm_control,tmp4$H20)
}

treatments<-c('DEX','TUNIC')

for (i in 1:length(treatments)){
	tmp<-lfsr[,which(names(lfsr)==treatments[i] | names(lfsr)=='ETOH' | names(lfsr)=='id')]
	tmp2<-subset(tmp,tmp[,1]<0.1 | tmp[,2]<0.1)
	tmp3<-pm[tmp2$id,which(names(pm)==treatments[i] | names(pm)=='ETOH' | names(pm)=='id' )]
	tmp3$diff<-log2(tmp3[,1]/tmp3[,2])
	tmp4<-subset(tmp3,diff < -1 | diff >1)
	out_diff<-c(out_diff,tmp4$diff)
	out_id<-c(out_id,tmp4$id)
	out_treat<-c(out_treat,as.character(rep(treatments[i],dim(tmp4)[1])))
	out_pm_treat<-c(out_pm_treat,tmp4[,which(names(tmp4)==treatments[i])])
	out_pm_control<-c(out_pm_control,tmp4$ETOH)
	}

out<-as.data.frame(cbind(out_id,out_diff,out_pm_treat,out_pm_control))
out$out_treat<-out_treat
out2<-merge(out,lfsr[,-c(1:12)],by.y='id',by.x='out_id')
out2$stronger_treatment<-abs(out2$out_pm_treat)>abs(out2$out_pm_control)
write.table(as.data.frame(out2),'10Jun21_eQTL_SNPs_response.txt',row.names=F,sep='\t')

as.matrix(table(out$out_treat))
tmp<-unique(out2[,c('out_treat','gene')])
as.matrix(table(tmp$out_treat))
