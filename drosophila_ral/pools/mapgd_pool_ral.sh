#!/bin/bash

module load SAMtools/1.6-foss-2016b
module load MAPGD/0.4.26

path=$1 #path to input
name=$2 #name of pool
alpha=$3 #value for each of the three different thresholds (22, 11 and 6) for the log-likelihood ratio (LLR) test (-a flag)

#sort and Index
samtools sort $path/$name_final.bam > $path/$name_final.sorted.bam
samtools index $path/$name_final.sorted.bam  

#run MAPGD
samtools view -H $path/$name_final.sorted.bam > $path/$name.header
mapgd proview -i $path/$name_final.mpileup -H $path/$name.header -o $path/$name
mapgd pool -i $path/$name.pro -a $alpha -o $path/$name.$alpha.allelefrequency-filtered