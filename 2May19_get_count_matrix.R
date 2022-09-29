#!/bin/bash Rscript

# NOTE: "FILEINFO" was replaced with information about each particular file

library(data.table)
ad=fread('/Genomics/ayroleslab2/alea/longevity/NovaSeq_Apr19/count_info/joint_call_FILEINFO.AD.filt.txt')
dp=fread('/Genomics/ayroleslab2/alea/longevity/NovaSeq_Apr19/count_info/joint_call_FILEINFO.DP.filt.txt')
# not sure why there are some duplicate entries?
ad=unique(ad)
dp=unique(dp)

ad$loc_id<-paste(ad$V1,ad$V2,sep='_')
dp$loc_id<-paste(dp$V1,dp$V2,sep='_')

both=merge(ad,dp,by='loc_id')
both2=subset(both,(V3.x+V4)==V3.y)

ad2=ad[which(ad$loc_id %in% both2$loc_id),]

setkey(ad2,V1,V2)
names=unique(ad2[,V2])
sites=unique(ad2[,V1])

counts=as.data.table(unique(ad2[,V1]))
setnames(counts,"V1","site")
for(i in 1:length(names))
{
  counts[,names[i] := ad2[.(as.factor(counts[,site]),names[i]),V3]]
#print(i) 
}
ref_counts<-counts

counts=as.data.table(unique(ad2[,V1]))
setnames(counts,"V1","site")
for(i in 1:length(names))
{
  counts[,names[i] := ad2[.(as.factor(counts[,site]),names[i]),V4]]
#print(i)
}
alt_counts<-counts

#both<-alt_counts[,-1]+ref_counts[,-1]
#missing_per_fly<-apply(both,2,function(x) length(which(is.na(x))))
#missing_per_site<-apply(both,1,function(x) length(which(is.na(x))))

#remove1<-(which(missing_per_fly/dim(both)[1]>0.9))
#remove2<-(which(missing_per_site/dim(both)[2]>0.75))

write.table(alt_counts,'/Genomics/ayroleslab2/shared/longevity/alternate_counts_FILEINFO.txt',row.names=F,sep='\t',quote=F)
write.table(ref_counts,'/Genomics/ayroleslab2/shared/longevity/reference_counts_FILEINFO.txt',row.names=F,sep='\t',quote=F)

