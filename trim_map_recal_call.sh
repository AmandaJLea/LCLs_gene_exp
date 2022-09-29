#!/bin/sh

# NOTE: "SAMPLEINFO" was replaced with the sample name

#######
# set parameters and directories
#######

module load samtools
module load java

cpus=4
ref_DIR=/Genomics/ayroleslab2/alea/ref_genomes/dmel/dmel-all-chromosome-r6.23.fasta
fastq_DIR=/Genomics/ayroleslab2/alea/archive_raw_fastq/Project_AYR_13970_B01_NAN_Lane.2019-04-19/Sample_SAMPLEINFO/fastq
picard_DIR=/Genomics/grid/users/alea/programs/picard-tools-1.141/picard.jar
sambamba_DIR=/Genomics/grid/users/alea/programs/sambamba_v0.6.6
GATK_DIR=/Genomics/grid/users/alea/programs/GenomeAnalysisTK.jar
variants_DIR=/Genomics/ayroleslab2/alea/longevity/NovaSeq_Dec18/joint_vcfs/joint_call_n2592.allchr.PASS.bed

tmp_DIR=/scratch/tmp/alea/longevity

R1=$fastq_DIR/SAMPLEINFO_*001.R1.fastq.gz
R2=$fastq_DIR/SAMPLEINFO_*001.R2.fastq.gz

R1_trim=$tmp_DIR/SAMPLEINFO.trim.R1.fastq.gz
R2_trim=$tmp_DIR/SAMPLEINFO.trim.R2.fastq.gz

#######
# trim
#######
cutadapt -e 0.1 --overlap 2 -a AGATCGGAAGAG -A AGATCGGAAGAG --minimum-length=20 --trim-n -o $R1_trim -p $R2_trim $R1 $R2

#######
# map
#######
bwa mem -M -t $cpus $ref_DIR $R1_trim $R2_trim | samtools view -Sbq 1 | samtools sort -@ $cpus -m 3G > $tmp_DIR/SAMPLEINFO_sort_uniq.bam

rm $R1_trim
rm $R2_trim

#######
# dedup
#######
java -Xmx40g -jar $picard_DIR MarkDuplicates I=$tmp_DIR/SAMPLEINFO_sort_uniq.bam O=$tmp_DIR/SAMPLEINFO_sort_uniq_dedup.bam M=$tmp_DIR/SAMPLEINFO_sort_uniq_dedup.txt

$sambamba_DIR index -t $cpus $tmp_DIR/SAMPLEINFO_sort_uniq_dedup.bam

#######
# add read groups
#######
java -jar $picard_DIR AddOrReplaceReadGroups I=$tmp_DIR/SAMPLEINFO_sort_uniq_dedup.bam O=$tmp_DIR/SAMPLEINFO_sort_uniq_dedup_rg.bam SO=coordinate RGLB=SAMPLEINFO RGPL=illumina RGPU=longevity RGSM=SAMPLEINFO 

$sambamba_DIR index -t $cpus $tmp_DIR/SAMPLEINFO_sort_uniq_dedup_rg.bam

#######
# base recalibration and application
#######
java -Xmx40g -jar $GATK_DIR -nct $cpus -R $ref_DIR -I $tmp_DIR/SAMPLEINFO_sort_uniq_dedup_rg.bam -T BaseRecalibrator -o $tmp_DIR/SAMPLEINFO.recal.table -knownSites $variants_DIR

java -Xmx40g -jar $GATK_DIR -nct $cpus -R $ref_DIR -T PrintReads -I $tmp_DIR/SAMPLEINFO_sort_uniq_dedup_rg.bam -BQSR $tmp_DIR/SAMPLEINFO.recal.table -o $tmp_DIR/SAMPLEINFO_recal.bam

$sambamba_DIR index -t $cpus $tmp_DIR/SAMPLEINFO_recal.bam

#######
# individual variant calling
#######
java -Xmx40g -jar $GATK_DIR -T HaplotypeCaller -nct $cpus -R $ref_DIR -I $tmp_DIR/SAMPLEINFO_recal.bam -ERC GVCF -o $tmp_DIR/SAMPLEINFO.recal.vcf -variant_index_type LINEAR -variant_index_parameter 128000 -maxAltAlleles 2

bgzip $tmp_DIR/SAMPLEINFO.recal.vcf
tabix -p vcf $tmp_DIR/SAMPLEINFO.recal.vcf.gz

#######
# clean up
#######
rm $tmp_DIR/SAMPLEINFO_sort_uniq*


