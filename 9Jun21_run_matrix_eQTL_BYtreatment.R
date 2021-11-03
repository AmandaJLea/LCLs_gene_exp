# for f in `cat chrom_names.txt` ; do cat matrix_eqtl.R | sed -e s/CHROMNAME/$f/g > matrix_eqtl.$f.R; done
# for f in `cat chrom_names.txt` ; do cat matrix_eqtl.sh | sed -e s/CHROMNAME/$f/g > matrix_eqtl.$f.sh; done
# rm commands1.sh; touch commands1.sh
# for f in `cat chrom_names.txt` ; do echo "sh "matrix_eqtl.$f.sh >>commands1.sh; done 
# sbatch -a 1-22%22 array1.sh

library(MatrixEQTL)
library(data.table)

# read in files
setwd("/Genomics/ayroleslab2/alea/LCL_stim/final_dataset")

names<-read.delim('igsr-1000 genomes 30x on grch38 - focal IDs.tsv',stringsAsFactors=F)
genepos <- read.delim("/Genomics/ayroleslab2/alea/ref_genomes/hg38/hg38_PCgene_start_end.txt", header = F)
genes=read.delim('31Mar21_all_runs_voom_norm_geneIDs.txt')
genoFile=fread('/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_20201028_3202_raw_GT/CCDG_14151_B01_GRM_WGS_2020-08-05_CHROMNAME.filtered.shapeit2.LCLs.traw')
exprFile=fread('31Mar21_all_runs_voom_resid.txt')
info=read.delim('31Mar21_filtered_sample_info_v2.txt',stringsAsFactors=F)
info$exp<-names(exprFile)
pca=read.delim('/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_20201028_3202_raw_GT/CCDG_14151_B01_GRM_WGS_2020-08-05_CHROMNAME.filtered.shapeit2.LCLs.eigenvec',header=F)
pca$id<-paste(pca$V1,pca$V2,sep='_')

# subset to treatment and get geno names for each sample
names1<-names[which(names$Sample.name %in% info$Individual.ID),]
names2<-names[which(names$Sample.GMme2 %in% info$Individual.ID),]
names3<-unique(rbind(names,names))

for (z in c('ACRYL','BAFF','BPA','ETOH','FSL1','GARD','IFNG','IGF','PFOA','DEX','TUNIC','H20') ){

info2<-subset(info, ( treatment==z) & (info$Individual.ID %in% names3$Sample.GMme2 | info$Individual.ID %in% names3$Sample.name))
info2$geno_ID<-NA

for (i in 1:dim(info2)[1]){
info2$geno_ID[i]<-as.character(ifelse(info2$Individual.ID[i] %in% names3$Sample.GMme2, names3$Sample.name[which(names3$Sample.GMme2 == info2$Individual.ID[i])] ,NA)) }

for (i in which(is.na(info2$geno_ID))){
info2$geno_ID[i]<-as.character(ifelse(info2$Individual.ID[i] %in% names3$Sample.name, names3$Sample.name[which(names3$Sample.name == info2$Individual.ID[i])] ,NA)) }

# file rearranging - geno
info2$geno_ID<-paste(info2$geno_ID,info2$geno_ID,sep="_")
genoFile2=as.data.frame(genoFile)
tmp<-subset(as.data.frame(table(genoFile2$SNP)),Freq==1)
genoFile3=as.data.frame(genoFile2[which(genoFile2$SNP %in% tmp$Var1),info2$geno_ID])
rownames(genoFile3)<-genoFile2$SNP[which(genoFile2$SNP %in% tmp$Var1)]

# file rearranging - SNP pos
library(stringr)
tmp2<-as.data.frame(str_split_fixed(genoFile2$SNP[which(genoFile2$SNP %in% tmp$Var1)], ":", 2))
tmp2$chr<-paste('chr',tmp2$V1,sep='')
SNPpos<-as.data.frame(cbind(rownames(genoFile3),tmp2[,c('chr')]))
SNPpos$loc<-tmp2$V2
write.table(SNPpos,'temp_CHROMNAME.txt',row.names=F,sep='\t')
SNPpos=read.delim('temp_CHROMNAME.txt')

# make covariate file
info3<-merge(info2,pca,by.x='geno_ID',by.y='id')
info3<-info3[order(info3$id),]

covars <- (info3[,c('V3','V4','V5','V6','V7')])
cvrt = SlicedData$new();
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$CreateFromMatrix(t(covars))

# make geno file
genos = SlicedData$new();
genos$fileOmitCharacters = "NA"; # denote missing values;
genos$fileSkipRows = 1;          # one row of column labels
genos$fileSkipColumns = 1;       # one column of row labels
genos$CreateFromMatrix(as.matrix(genoFile3))

# make gene loc file
genepos2=genepos[,c('V4','V5','V1','V2','V3')]
genes$V2<-1:dim(genes)[1]
genepos2=merge(genepos2,genes,by.y='V1',by.x='V4')
genepos2=genepos2[order(genepos2$V2.y),]
genepos2=genepos2[!duplicated(genepos2$V5), ]

# make exp file
exprFile=as.data.frame(exprFile)
fc=exprFile[genepos2$V2.y,(info2$exp)]
rownames(fc)<-genepos2$V5
genepos2=genepos2[,2:5]

expr = SlicedData$new();
expr$fileOmitCharacters = "NA"; # denote missing values;
expr$fileSkipRows = 0;          # one row of column labels
expr$fileSkipColumns = 0;       # one column of row labels
expr$CreateFromMatrix(as.matrix(fc))

# checks
identical(names(fc),info2$exp)
identical(as.character(genepos2$V5),rownames(fc))

########
# RUN
#########

i=0

 outFile_covar <- paste("31Mar21_LCLs_trans_eQTL_CHROMNAME.",z,i,".matEQTL",sep='')
 outFile_noCovar.cis <- paste("31Mar21_LCLs_cis_eQTL_CHROMNAME.",z,i,".matEQTL",sep='')
out_withCovar = Matrix_eQTL_main(
  cvrt = cvrt,
  snps = genos,
  gene = expr,
  output_file_name = outFile_covar,
 output_file_name.cis = outFile_noCovar.cis,
  snpspos = SNPpos,
  genepos = genepos2,
  pvOutputThreshold = 0,
  pvOutputThreshold.cis = 1,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  cisDist = 5e5,
  verbose = T,
  pvalue.hist = F,
  min.pv.by.genesnp = F,
  noFDRsaveMemory = F)
}
