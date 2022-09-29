#!/bin/sh

###########
# joint call by chromosome
###########

cpus=4
ref_DIR=/Genomics/ayroleslab2/alea/ref_genomes/dmel/dmel-all-chromosome-r6.23.fasta
out_DIR=/Genomics/ayroleslab2/alea/longevity/NovaSeq_Apr19/
GATK_DIR=/Genomics/grid/users/alea/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
variants_DIR=/Genomics/ayroleslab2/alea/longevity/NovaSeq_Dec18/joint_vcfs/joint_call_n2592.3R.PASS.vcf.gz

java -Xmx80g -Djava.io.tmpdir=$out_DIR -jar $GATK_DIR \
   -T GenotypeGVCFs \
   -R $ref_DIR -nt $cpus \
   --dbsnp $variants_DIR -V:VCF /Genomics/ayroleslab2/alea/longevity/NovaSeq_Apr19/scripts/vcfs.FILEINFO.list \
   -L 3R -o $out_DIR/joint_call_3R_FILEINFO.vcf --max_alternate_alleles 2

###########
# get filtered SNPs
###########

# https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set

java -Xmx80g -jar $GATK_DIR -R $ref_DIR -T SelectVariants --variant $out_DIR/joint_call_3R_FILEINFO.vcf -selectType SNP -o $out_DIR/joint_call_3R_FILEINFO.SNPs.vcf

rm $out_DIR/joint_call_3R_FILEINFO.vcf*

java -Xmx80g -Djava.io.tmpdir=/scratch/tmp/alea/longevity -jar $GATK_DIR -R $ref_DIR -T VariantFiltration -o $out_DIR/joint_call_3R_FILEINFO.SNPs.filt.vcf --variant $out_DIR/joint_call_3R_FILEINFO.SNPs.vcf --filterExpression "QUAL < 20.0" --filterExpression "QD < 2.0" --filterExpression "FS > 60.0" --filterExpression "MQ < 35.0" --filterExpression "MQRankSum < -12.5" --filterExpression "ReadPosRankSum < -8.0" --filterName "LowQual" --filterName "LowQD" --filterName "HighFS" --filterName "LowMQ" --filterName "LowMQRankSum" --filterName "LowReadPosRankSum"

rm $out_DIR/joint_call_3R_FILEINFO.SNPs.vcf*

~/programs/plink_1.90 --vcf $out_DIR/joint_call_3R_FILEINFO.SNPs.filt.vcf --allow-extra-chr --maf 0.01 --set-missing-var-ids @:#:\$1,\$2 --indep-pairwise 200 25 0.5 --geno 0.75 --biallelic-only strict --recode vcf --out $out_DIR/joint_call_3R_FILEINFO.SNPs.PASS

sed -e s/:/'\t'/g $out_DIR/joint_call_3R_FILEINFO.SNPs.PASS.prune.in | awk '{OFS="\t";print $1,$2-1,$2}' > $out_DIR/joint_call_3R_FILEINFO.SNPs.PASS.bed

bgzip $out_DIR/joint_call_3R_FILEINFO.SNPs.PASS.vcf
tabix -p vcf $out_DIR/joint_call_3R_FILEINFO.SNPs.PASS.vcf.gz

bcftools filter -e "QUAL < 20" --threads 4 -R $out_DIR/joint_call_3R_FILEINFO.SNPs.PASS.bed -O z -o $out_DIR/joint_call_3R_FILEINFO.SNPs.PASS2.vcf.gz $out_DIR/joint_call_3R_FILEINFO.SNPs.PASS.vcf.gz

rm $out_DIR/joint_call_3R_FILEINFO.SNPs.PASS.*

# keep some unfiltered SNPs and final filtered set
bgzip $out_DIR/joint_call_3R_FILEINFO.SNPs.filt.vcf
tabix -p vcf $out_DIR/joint_call_3R_FILEINFO.SNPs.filt.vcf.gz

tabix -p vcf $out_DIR/joint_call_3R_FILEINFO.SNPs.PASS2.vcf.gz


