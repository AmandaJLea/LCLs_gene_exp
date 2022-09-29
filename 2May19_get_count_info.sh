#!/bin/bash

# NOTE: "FILEINFO" was replaced with information about each particular file

out_DIR=/Genomics/ayroleslab2/alea/longevity/NovaSeq_Apr19/joint_vcfs
vcf_file=$out_DIR/joint_call_FILEINFO.SNPs.filt.vcf.gz

out_file1=/Genomics/ayroleslab2/alea/longevity/NovaSeq_Apr19/count_info/joint_call_FILEINFO.DP.txt
out_file2=/Genomics/ayroleslab2/alea/longevity/NovaSeq_Apr19/count_info/joint_call_FILEINFO.AD.txt

bcftools query -f '[%CHROM:%POS\t%SAMPLE\t%DP\n]' --regions-file $out_DIR/joint_call_FILEINFO.SNPs.PASS2.vcf.gz $vcf_file | sed 's/,/\t/g'> $out_file1
bcftools query -f '[%CHROM:%POS\t%SAMPLE\t%AD\n]' --regions-file $out_DIR/joint_call_FILEINFO.SNPs.PASS2.vcf.gz $vcf_file | sed 's/,/\t/g'> $out_file2

# remove sites with no counts
awk '{OFS="\t"; if ($3>0) print $1,$2,$3}' $out_file1 > /Genomics/ayroleslab2/alea/longevity/NovaSeq_Apr19/count_info/joint_call_FILEINFO.DP.filt.txt

awk '{OFS="\t"; if ($3>0 || $4>0) print $1,$2,$3,$4}' $out_file2 > /Genomics/ayroleslab2/alea/longevity/NovaSeq_Apr19/count_info/joint_call_FILEINFO.AD.filt.txt

rm $out_file1
rm $out_file2

