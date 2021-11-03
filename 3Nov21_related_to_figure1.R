################
# map of study samples
################

setwd("~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL")

gps=read.delim('igsr_populations_wGPS.tsv')
gps2<-subset(gps,Include==1)

tmp<-as.data.frame(unique(gps2[,c('Superpopulation.code','Population.code')]))
tmp<-tmp[order(tmp$Superpopulation.code,tmp$Population.code),]
gps2$Population <- factor(gps2$Population.code, levels = tmp$Population.code)

library(maps);library(mapdata)
library(ggplot2);library(ggmap); library(devtools)
library(ggrepel)

# sampling locations - zoomed in
register_google('AIzaSyBqRWaHZTtsvaSOIwDdHjmTVPD-cXvaKrE')

map2 <- get_map(location = c(lon = mean(gps2$Population.longitude), lat = mean(gps2$Population.latitude)),maptype = "watercolor",zoom=2,color="bw")
m2 <- ggmap(map2)

m2 +  geom_point(data=gps2,aes(x=Population.longitude, y=Population.latitude,colour=Population), size=4) + geom_point(data=gps2,aes(x=Population.longitude, y=Population.latitude),shape = 1,size = 4,colour = "black")+ scale_colour_brewer(palette = "Paired")+theme_bw(15)+xlab("Longitude")+ylab("Latitude")+geom_label_repel(data=gps2,aes(x=Population.longitude, y=Population.latitude,label=Population),box.padding   = 0.35, point.padding = 0.5)
    
################
# PCA of genotype data
################

library(data.table)

# read in files
setwd("/Genomics/ayroleslab2/alea/LCL_stim/final_dataset")

names<-read.delim('igsr-1000 genomes 30x on grch38 - focal IDs.tsv',stringsAsFactors=F)
pca=read.delim('/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_20201028_3202_raw_GT/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2.LCLs.eigenvec',header=F)

pca2<-merge(names,pca,by.x='Sample.name',by.y='V1')
write.table(pca2,"10Sep21_genotype_PCs.txt",row.names=F,sep='\t')

####

setwd("~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL")
data=read.delim('10Sep21_genotype_PCs.txt')

tmp<-as.data.frame(unique(data[,c('Superpopulation.code','Population.code')]))
tmp<-tmp[order(tmp$Superpopulation.code,tmp$Population.code),]
data$Population <- factor(data$Population.code, levels = tmp$Population.code)

library(ggplot2)
ggplot(data, aes(x=V3, y=V4,col=Population)) + geom_point()+scale_colour_brewer(palette = "Paired")+theme_bw(15)+xlab('PC 1')+ylab('PC 2')

################
# treatment effects
################

setwd("~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL")
library(mashr)
library(gplots)
library(corrplot)

lfsr=read.delim('31Mar21_LFSR_DE_mashr.txt')
pm=read.delim('31Mar21_pm_DE_mashr.txt')

names(pm)[5]<-'FSL-1'
names(pm)[8]<-'IGF-1'
corrplot(cor(pm),method='ellipse',type='upper')  

tmp<-apply(lfsr,2,function(x) length(which(x<0.05)))
# assess sharing
lfsr$sig<-apply(lfsr,1,function(x) length(which(x<0.1)))
lfsr$id<-1:dim(lfsr)[1]
lfsr_sig1<-subset(lfsr,sig>0 )

lfsr_sig1$within2_v2<-NA
lfsr_sig1$focal<-NA
lfsr_sig1$info<-NA

within2_info1<-c();within2_info2<-c()

for (i in 1:dim(lfsr_sig1)[1]){
	z<-min(lfsr_sig1[i,1:11])[1]
	lfsr_sig1$focal[i]<-which(lfsr_sig1[i,1:11] == z)
	focal<-pm[lfsr_sig1$id[i],which(lfsr_sig1[i,1:11] == z)][1]
	not_focal<-pm[lfsr_sig1$id[i],which(lfsr_sig1[i,1:11] != z)]
	not_focal2<-log2(not_focal/as.numeric(focal))
	lfsr_sig1$within2_v2[i]<- length(which(not_focal2 > -1 & not_focal2<1))
	within2_info1<- c(within2_info1,c(which(lfsr_sig1[i,1:11] == z)[1],(which(not_focal2 > -1 & not_focal2<1))))
	within2_info2<- c(within2_info2,rep(lfsr_sig1$id[i],1+length(which(not_focal2 > -1 & not_focal2<1))))
}

# what's going on with 1-2 categories?
tmp<-subset(lfsr_sig1,within2_v2==0)
table(tmp$focal)
tmp<-subset(lfsr_sig1,within2_v2==1)
table(tmp$focal)
tmp<-subset(lfsr_sig1,within2_v2==2)
table(tmp$focal)

tmp2<-as.data.frame(cbind(within2_info1,within2_info2))
tmp3<-subset(tmp2,within2_info2 %in% tmp$id)
x<-as.data.frame(table(tmp3[,1]))
ggplot(x, aes(x=Var1,y=Freq )) +geom_bar(stat = "identity",fill="steelblue" )+theme_bw(13)

# remove DEX
lfsr_noDEX<-lfsr[,-10]
# assess sharing
lfsr_noDEX$sig<-apply(lfsr_noDEX[,1:10],1,function(x) length(which(x<0.1)))
lfsr_noDEX$id<-1:dim(lfsr_noDEX)[1]
lfsr_sig1_noDEX<-subset(lfsr_noDEX,sig>0 )

pm_noDEX<-pm[,-10]
lfsr_sig1_noDEX$within2_v2<-NA
lfsr_sig1_noDEX$focal<-NA

for (i in 1:dim(lfsr_sig1_noDEX)[1]){
	z<-min(lfsr_sig1_noDEX[i,1:10])[1]
	lfsr_sig1_noDEX$focal[i]<-which(lfsr_sig1_noDEX[i,1:10] == z)
	focal<-pm_noDEX[lfsr_sig1_noDEX$id[i],which(lfsr_sig1_noDEX[i,1:10] == z)][1]
	not_focal<-pm_noDEX[lfsr_sig1_noDEX$id[i],which(lfsr_sig1_noDEX[i,1:10] != z)]
	not_focal2<-log2(not_focal/as.numeric(focal))
	lfsr_sig1_noDEX$within2_v2[i]<- length(which(not_focal2 > -1 & not_focal2<1))
}

x<-as.data.frame(table(lfsr_sig1_noDEX$within2_v2+1))
ggplot(x, aes(x=Var1,y=Freq )) +geom_bar(stat = "identity",fill="steelblue" )+theme_bw(13)

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

heatmap.2(as.matrix(out),dendrogram = "none",trace="none", notecol = "Black",notecex = 2,key=T,col=colorRampPalette(c("Steel Blue","White","Red")),margins=c(12,9),density.info='none')

apply(lfsr,2,function(x) length(which(x<0.1)))
