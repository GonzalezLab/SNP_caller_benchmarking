#!/bin/bash  

module load Python/2.7.12-foss-2016b 
module load SAMtools/1.6-foss-2016b

path=$1 #path to input
name=$2 #name of pool

#run PoolSNP
$path/PoolSNP/PoolSNP.sh \
mpileup=$path$name_final.mpileup \
output=$path/$name.output \
reference=$path/dmel_6.12.fa \
names=$name \
min-cov=8 \
max-cov=0.95 \
min-count=2 \
min-freq=0.016 \
bq=15 \
jobs=1 \
BS=1
