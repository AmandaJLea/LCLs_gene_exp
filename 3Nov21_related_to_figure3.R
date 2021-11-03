############
# change in h2
############

setwd("~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL")

treat1<-c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA','DEX','TUNIC','H20')

out=read.delim('25Oct21_GCTA_ACRYL_results_1GRM.txt',header=F)
out$info<-'ACRYL'
for (k in 2:12){
tmp=read.delim(paste('25Oct21_GCTA_',treat1[k],'_results_1GRM.txt',sep=''),header=F)
tmp$info<-treat1[k]
out=rbind(out,tmp)
}
out_subsamp<-out

out=read.delim('GCTA_ACRYL_results_1GRM.txt',header=F)
out$info<-'ACRYL'
for (k in 2:12){
tmp=read.delim(paste('GCTA_',treat1[k],'_results_1GRM.txt',sep=''),header=F)
tmp$info<-treat1[k]
out=rbind(out,tmp)
}

treat1_df<-as.data.frame(treat1)
treat1_df$treat2<-treatment<-c(1,1,1,0,1,1,1,1,1,1,1,0)
out2=merge(treat1_df,out_subsamp,by.x='treat1',by.y='info')

tmp<-as.data.frame(aggregate(out$V2~out$info,FUN=mean))
names(tmp)<-c('condition','mean')
tmp$se<-as.data.frame(aggregate(out$V2~out$info,FUN=function(x) sd(x)/sqrt(10156)))[,2]
tmp$condition2<-c("ACRYL","BAFF",  "BPA" ,  "DEX" ,  "ETOH" , "FSL-1",  "GARD",  "H20" ,  "IFNG" ,"IGF-1" ,  "PFOA",  "TUNIC")
tmp$treatment<-c(1,1,1,1 ,  0, 1,1,0  ,1,1,1,1)

tmp2<-as.data.frame(aggregate(out_subsamp$V2~out_subsamp$info,FUN=mean))
names(tmp2)<-c('condition','mean')
tmp2$se<-as.data.frame(aggregate(out_subsamp$V2~out_subsamp$info,FUN=function(x) sd(x)/sqrt(10156)))[,2]
tmp2$condition2<-c("ACRYL","BAFF",  "BPA" ,  "DEX" ,  "ETOH" , "FSL-1",  "GARD",  "H20" ,  "IFNG" ,"IGF-1" ,  "PFOA",  "TUNIC")
tmp2$treatment<-c(1,1,1,1 ,  0, 1,1,0  ,1,1,1,1)

ggplot(tmp, aes(x=factor(treatment), y=mean,color=condition2)) +
   geom_point(position=position_dodge(), stat="identity") +
geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)+theme_bw(13)+scale_color_brewer(palette="Paired")+ylim(0.05,0.4)

ggplot(tmp2, aes(x=factor(treatment), y=mean,color=condition2)) +
   geom_point(position=position_dodge(), stat="identity") +
geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)+theme_bw(13)+scale_color_brewer(palette="Paired")+ylim(0.05,0.4)

ggplot(tmp, aes(x=factor(treatment), y=mean)) + geom_violin()+geom_boxplot(outlier.shape = NA,width=0.2) + theme_bw(13)+ scale_color_brewer(palette="Paired")+theme(legend.position = "top")+geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.3), aes(color=factor(condition2)))+ylim(0.05,0.4)

ggplot(tmp2, aes(x=factor(treatment), y=mean)) + geom_violin()+geom_boxplot(outlier.shape = NA,width=0.2) + theme_bw(13)+ scale_color_brewer(palette="Paired")+theme(legend.position = "top")+geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.3), aes(color=factor(condition2)))+ylim(0.05,0.4)

############
# Number of eQTL in each condition and sharing
############

setwd("~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL")
library(mashr)
library(gplots)
library(corrplot)

lfsr=read.delim('31Mar21_LFSR_eQTL_mashr_nocor.txt')
pm=read.delim('31Mar21_pm_eQTL_mashr_nocor.txt')
sharing=read.delim('9Jun21_analyses/10Jun21_eQTL_SNPs_sharing.txt')
    
x<-as.data.frame(table(sharing$within2_v2+1))
ggplot(x, aes(x=Var1,y=Freq )) +geom_bar(stat = "identity",fill="steelblue" )+theme_bw(13)
x2<-as.data.frame(table(unique(sharing[,c('within2_v2','gene')])$within2_v2+1))
ggplot(x2, aes(x=Var1,y=Freq )) +geom_bar(stat = "identity",fill="steelblue" )+theme_bw(13)

lfsr$sig<-apply(lfsr,1,function(x) length(which(x<0.1)))
lfsr$id<-1:dim(lfsr)[1]
lfsr_sig1<-subset(lfsr,sig>0 )

# SNP level sharing
out<-as.data.frame(matrix(ncol=12,nrow=12))
rownames(out)<-colnames(out)<-names(lfsr_sig1)[1:12]

for (i in c(names(lfsr_sig1)[1:12])){
	for (k in c(names(lfsr_sig1)[1:12])){
	lfsr_tmp<-lfsr[,c(k,i)]
	pm_tmp<-pm[,c(k,i)]
	lfsr_tmp$sig<-apply(lfsr_tmp,1,function(x) length(which(x<0.1)))
	lfsr_tmp$id<-1:dim(lfsr_tmp)[1]
	lfsr_tmp2<-subset(lfsr_tmp,sig>0)
	pm_tmp$id<-1:dim(pm_tmp)[1]
	both<-merge(lfsr_tmp2,pm_tmp,by='id')
	both$log2<-log2(both[,5]/both[,6])
	out[i,k]<-dim(subset(both,log2> -1 & log2<1))[1] }}

rownames(out)[6]<-colnames(out)[6]<-'FSL-1'
rownames(out)[10]<-colnames(out)[10]<-'IGF-1'

heatmap.2(as.matrix(out),dendrogram = "none",trace="none", notecol = "Black",notecex = 2,key=T,col=colorRampPalette(c("Steelblue","White","Red")),margins=c(12,9),density.info='none')

