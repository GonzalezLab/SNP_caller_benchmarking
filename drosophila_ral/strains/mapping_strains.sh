#!/bin/bash 

path=$1 #full path to reference genome 
read1=$2 #full path to read 1 input
read2=$3 #full path to read 2 input
name=$4 #name of the output file
out=$5 #full path to the output 

module load cutadapt/1.9.1-foss-2016b-Python-2.7.12
module load BWA/0.7.16a-foss-2016b
module load SAMtools/1.6-foss-2016b
module load GATK/3.4-46
module load GATK/3.7-Java-1.8.0_74
module load picard/2.8.3

##trim sequences
cutadapt \
-q 18 \
-o $out/trimmed-$name-1.fq.gz \
-p $out/trimmed-$name-2.fq.gz \
$read1 \
$read2

mkdir $out/mapping_bwa

#Built an index
bwa index -a is $path/dmel_6.12.fa

##mapping
bwa mem -M -t 5 $path/dmel_6.12.fa \
	$out/trimmed-$name-1.fq.gz $out/trimmed-$name-2.fq.gz | samtools view -Sbh -q 20 -F 0x100 - > $out/mapping_bwa/$name.bam

mkdir $out/mapping_bwa/dedup-report

##sort
java -Xmx20G -XX:+UseSerialGC -Dsnappy.disable=true -jar $PICARDPATH/picard.jar \
	SortSam I=$out/mapping_bwa/$name.bam O=$out/mapping_bwa/$name-sort.bam SO=coordinate VALIDATION_STRINGENCY=SILENT

##de-duplicate
java -Xmx20g -XX:+UseSerialGC -Dsnappy.disable=true -jar $PICARDPATH/picard.jar \
	MarkDuplicates REMOVE_DUPLICATES=true I=$out/mapping_bwa/$name-sort.bam O=$out/mapping_bwa/$name-dedup.bam M=$out/mapping_bwa/dedup-report/$name.txt VALIDATION_STRINGENCY=SILENT

##create DICT file
java -Xmx20G -XX:+UseSerialGC -jar $PICARDPATH/picard.jar CreateSequenceDictionary \
	REFERENCE=$path/dmel_6.12.fa \
	OUTPUT=$path/dmel_6.12.dict

##create FASTA index file
samtools faidx $path/dmel_6.12.fa

##add read group to BAM files
java -Xmx10G -XX:+UseSerialGC -jar $PICARDPATH/picard.jar AddOrReplaceReadGroups \
	INPUT=$out/mapping_bwa/$name-dedup.bam OUTPUT=$out/mapping_bwa/$name-dedup_rg.bam \
	SORT_ORDER=coordinate RGID=$name RGLB=$name RGPL=illumina RGSM=sample RGPU=name \
	CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT

##remove BAM without ReadGroup
rm $out/mapping_bwa/$name-dedup.bam

##generate target List of InDel positions 
mkdir $out/mapping_bwa/realign_list
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $path/dmel_6.12.fa \
	-I $out/mapping_bwa/$name-dedup_rg.bam -o $out/mapping_bwa/realign_list/$name.list

##re-align around InDels
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $path/dmel_6.12.fa \
	-I $out/mapping_bwa/$name-dedup_rg.bam -targetIntervals $out/mapping_bwa/realign_list/$name.list \
	-o $out/mapping_bwa/$name-dedup_rg_InDel.bam

##filter out SNPs in sex chr and chr 4
samtools sort $name-dedup_rg_InDel.bam > $name_sorted.bam
samtools index $name_sorted.bam

samtools view $name_sorted.bam chr2L -b > $name_sorted_2L.bam
samtools view $name_sorted.bam chr2R -b > $name_sorted_2R.bam
samtools view $name_sorted.bam chr3L -b > $name_sorted_3L.bam
samtools view $name_sorted.bam chr3R -b > $name_sorted_3R.bam

samtools merge $name_final.bam $name_sorted_2L.bam $name_sorted_2R.bam $name_sorted_3L.bam $name_sorted_3R.bam
