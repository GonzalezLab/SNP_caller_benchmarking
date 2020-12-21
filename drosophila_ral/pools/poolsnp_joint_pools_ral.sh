#!/bin/bash  

module load Python/2.7.12-foss-2016b 
module load SAMtools/1.6-foss-2016b

path=$1 #path to input
name1=$2 #name of pool1
name2=$3 #name of pool2
mf=$4 #miss-fraction (0.10 or 0.80)

#Generate mpileup
samtools mpileup -f -B $path/dmel_6.12.fa \
	$path/$name1_final.bam \
	$path/$name2_final.bam \
    > $path/$name1.$name2.mpileup

#run PoolSNP
$path/PoolSNP/PoolSNP.sh \
mpileup=$path/$name1.$name2.mpileup \
output=$path/$name1.$name2.mf.$mf.output \
reference=$path/dmel_6.12.fa \
names=$name1,$name2 \
min-cov=8 \
max-cov=0.95 \
min-count=2 \
min-freq=0.016 \
miss-frac=$mf \
bq=15 \
jobs=1 \
BS=1
