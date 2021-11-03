#######
# GWAS enrichment - get data
#######

cd /Genomics/ayroleslab2/alea/LCL_stim/final_dataset/

library(data.table)
GWAS=fread('/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GWAS_catalog.txt')
GWAS=unique(GWAS[,c(8,12,13,29),with=F])
GWAS=subset(GWAS,PVALUE_MLOG>8)
GWAS$SNP<-paste(GWAS$CHR_ID,GWAS$CHR_POS,sep=':')

eQTL=fread('/Genomics/ayroleslab2/alea/LCL_stim/final_dataset/10Jun21_eQTL_SNPs_sharing.txt')
eQTL$SNP<-paste('chr',eQTL$SNP,sep='')
tested_SNPs=fread('31Mar21_LCLs_all_tested_SNPs.txt')
tested_SNPs2=as.data.frame(unique(tested_SNPs$SNP))
names(tested_SNPs2)<-'SNP'
tested_SNPs2$SNP2<-paste('chr',tested_SNPs2$SNP,sep='')

tested_SNPs2$iseQTL<-0
tested_SNPs2$iseQTL[which(tested_SNPs2$SNP2 %in% subset(eQTL,within2_v2> -1 )$SNP)]<-1

tested_SNPs2$isshared<-0
tested_SNPs2$isshared[which(tested_SNPs2$SNP2 %in% subset(eQTL,within2_v2==11 )$SNP)]<-1

tested_SNPs2$isGxE<-0
tested_SNPs2$isGxE[which(tested_SNPs2$SNP2 %in% subset(eQTL,within2_v2!=11 )$SNP)]<-1

tested_SNPs2$isGWAS<-0
tested_SNPs2$isGWAS[which(tested_SNPs2$SNP %in% GWAS$SNP)]<-1

write.table(tested_SNPs2,"13Sep21_GWAS_results.txt",row.names=F,sep='\t')

#######
# LOF enrichment - get data and results
#######

setwd('~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL')
library(data.table)

eQTL=fread('~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL/9Jun21_analyses/10Jun21_eQTL_SNPs_sharing.txt')
eQTL$SNP<-paste('chr',eQTL$SNP,sep='')
tested_SNPs=fread('31Mar21_LCLs_all_tested_SNPs.txt')
tested_SNPs2=as.data.frame(unique(tested_SNPs$SNP))
names(tested_SNPs2)<-'SNP'
tested_SNPs2$SNP2<-paste('chr',tested_SNPs2$SNP,sep='')
all_genes<-as.data.frame(table(tested_SNPs$gene))

lof_genes<-fread('nature19057-SI-Table-13.txt')
lof_genes<-lof_genes[,c('gene','pLI','clinvar_path','omim')]
lof_genes<-subset(lof_genes,pLI>0.9)

omim_genes<-read.delim('OMIM_genemap2.txt')

all_genes$is_lof<-0
all_genes$is_lof[which(all_genes$Var1 %in% lof_genes$gene)]<-1

all_genes$is_omim<-0
all_genes$is_omim[which(all_genes$Var1 %in% omim_genes$Gene.Symbols)]<-1

all_genes$is_shared<-0
all_genes$is_shared[which(all_genes$Var1 %in% subset(eQTL,within2_v2==11)$gene)]<-1

all_genes$is_GxE<-0
all_genes$is_GxE[which(all_genes$Var1 %in% subset(eQTL,within2_v2!=11)$gene)]<-1

all_genes$is_eGene<-0
all_genes$is_eGene[which(all_genes$Var1 %in% eQTL$gene)]<-1

# eQTL vs not eQTL
fisher.test(table(all_genes$is_lof,all_genes$is_eGene))
fisher.test(table(all_genes$is_omim,all_genes$is_eGene))

# remove overlapping GxE and eGenes
all_genes2<-subset(all_genes,(is_shared==1 & is_GxE==0) | (is_shared==0 & is_GxE==0))

# shared eGene vs not eGene
fisher.test(table(all_genes2$is_lof,all_genes2$is_shared))
fisher.test(table(all_genes2$is_omim,all_genes2$is_shared))

all_genes2<-subset(all_genes, (is_shared==0 & is_GxE==1) | (is_shared==0 & is_GxE==0))

# GxE eGene vs not eGene
fisher.test(table(all_genes2$is_lof,all_genes2$is_GxE))
fisher.test(table(all_genes2$is_omim,all_genes2$is_GxE))

######
# GWAS, iHS - get data and results
######

setwd("~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL")

library(data.table)
data=fread('13Sep21_GWAS_results.txt')

fisher.test(table(data$iseQTL,data$isGWAS))

tested_SNPs2=subset(data,(isshared==1 & isGxE==0 ) | (isshared==0 & isGxE==0 ) )
fisher.test(table(tested_SNPs2$isshared,tested_SNPs2$isGWAS))

tested_SNPs2=subset(data,(isshared==0 & isGxE==1 ) | (isshared==0 & isGxE==0 ) )
fisher.test(table(tested_SNPs2$isGxE,tested_SNPs2$isGWAS))

######
# PLOT enrichment
######

setwd('~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL/')
data=read.delim('14Sep21_enrichment_results.txt',header=T,stringsAsFactors=F)

data2=subset(data,is.na(Other))
data2$Category2 = factor(data2$Category, levels = c('eQTL','Ubiquitous eQTL','Context-dependent eQTL'))

library(ggplot2)

ggplot(data2, aes(x=Category2, y=Estimate, colour=Category2)) + 
    geom_errorbar(aes(ymin=Se1, ymax=Se2), width=.1) +geom_hline(yintercept=1,lty=2,col='grey')+
    geom_point()+facet_wrap(~Outcome)+scale_color_brewer(palette="Set1")+theme_bw(13)

data2=subset(data,Other!='NA')[1:11,]
data2$Other[which(data2$Other=='Novel contaminant or cell stressor ')]<-'Contaminant or stressor'
data2$Other2<-relevel(as.factor(data2$Other),ref='Immune stimulant')
summary(lm(data2$Estimate~data2$Other2))

ggplot(data2, aes(x=Category, y=Estimate)) + 
    geom_errorbar(aes(ymin=Se1, ymax=Se2), width=.1) +geom_hline(yintercept=1,lty=2,col='grey')+ geom_point()+facet_wrap(~Other,scales='free_x')+theme_bw(11)+ theme(axis.text.x = element_text(angle = 45,hjust=1))

######
# PLOT TWAS
######

setwd('~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL/9Jun21_analyses')

data_sig=read.delim('6Jul21_TWAS_eQTL_gene_results.txt',stringsAsFactors=F)
data_sig$info[which(data_sig$info=='shared')]<-'Ubiquitous'
data_sig$info[which(data_sig$info=='GxE')]<-'Context-dependent'

library(RColorBrewer)
brewer.pal(3, "Set1")

library(ggplot2)
ggplot(data_sig, aes(fill=info, y=(enrich), x=trait2)) + geom_bar(position="dodge", stat="identity")+theme_bw(13)+scale_fill_manual(values=c("#4DAF4A","#377EB8" ))+coord_flip()

######
# PLOT conservation
######

setwd('~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL')
library(data.table)

eQTL=fread('~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL/9Jun21_analyses/10Jun21_eQTL_SNPs_sharing.txt')
eQTL$SNP<-paste('chr',eQTL$SNP,sep='')
tested_SNPs=fread('31Mar21_LCLs_all_tested_SNPs.txt')
tested_SNPs2=as.data.frame(unique(tested_SNPs$SNP))
names(tested_SNPs2)<-'SNP'
tested_SNPs2$SNP2<-paste('chr',tested_SNPs2$SNP,sep='')
all_genes<-as.data.frame(table(tested_SNPs$gene))

all_genes$is_shared<-0
all_genes$is_shared[which(all_genes$Var1 %in% subset(eQTL,within2_v2==11)$gene)]<-1

all_genes$is_GxE<-0
all_genes$is_GxE[which(all_genes$Var1 %in% subset(eQTL,within2_v2!=11)$gene)]<-1

all_genes$is_eGene<-0
all_genes$is_eGene[which(all_genes$Var1 %in% eQTL$gene)]<-1

cons1=read.delim('ALLchr.phastCons100way.MeanByExon.bed')
cons2=read.delim('ALLchr.phyloP100way.MeanByExon.bed')
cons=merge(cons1,cons2,by='gene')

all_genes=merge(all_genes,cons,by.y='gene',by.x='Var1')

all_genes$cat<-NA

all_genes2<-subset(all_genes,(is_shared==0 & is_GxE==1) )
all_genes$cat[which(all_genes$Var1 %in% all_genes2$Var1)]<-'2_GxE_only'

all_genes2<-subset(all_genes,(is_shared==1 & is_GxE==0) )
all_genes$cat[which(all_genes$Var1 %in% all_genes2$Var1)]<-'1_shared_only'

all_genes2<-subset(all_genes,(is_shared==0 & is_GxE==0) )
all_genes$cat[which(all_genes$Var1 %in% all_genes2$Var1)]<-'0_not_eGene'

all_genes2<-subset(all_genes,(is_shared==1 & is_GxE==1) )
all_genes$cat[which(all_genes$Var1 %in% all_genes2$Var1)]<-'both_types_eGene'

tmp<-subset(all_genes,cat!='both_types_eGene')
library(ggplot2)
ggplot(tmp, aes(y= mean_score.x,x=factor(cat),fill=factor(cat))) + geom_violin() + geom_boxplot(width=0.1,outlier.shape = NA) +theme_bw(13)+ylab('Mean per-gene phastCons score')+ scale_fill_manual(values=c("Grey","#377EB8","#4DAF4A" ))+ylim(0,1)

ggplot(tmp, aes(y= mean_score.y,x=factor(cat),fill=factor(cat))) + geom_violin() + geom_boxplot(width=0.1,outlier.shape = NA) +theme_bw(13)+ylab('Mean per-gene phyloP score')+ scale_fill_manual(values=c("Grey","#377EB8","#4DAF4A" ))+ylim(0,5)



