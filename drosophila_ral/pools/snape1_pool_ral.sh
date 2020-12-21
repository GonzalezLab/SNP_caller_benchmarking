#!/bin/bash                                                                                                                                                                                     
module load snape-pooled/1 
module load R/3.5.1-foss-2018b

path=$1 #path to input
name=$2 #name of pool
nchr=$3 #number of chromosomes
 
#remove positions with 0 reads and run SNAPE with flat prior
if [ $(cut -f4 $path/$name_final.mpileup | grep -c "^0$") != 0  ]
then
	snape-pooled -nchr $nchr -theta 0.01 -D 0.1 -fold unfolded  -priortype flat < $path/$name_final.mpileup > $path/$name_final.$nchr.chr.snape.out
else
	cat $path/mapping_bwa/$name_final.mpileup | grep -v $'\t0\t' > $path/$name_final.mod.mpileup
	snape-pooled -nchr $nchr -theta 0.01 -D 0.1 -fold unfolded  -priortype flat < $path/$name_final.mod.mpileup > $path/snape/$name_final.$nchr.chr.snape.out
fi

#Filter SNPs with with a posterior probability < 0.90
cat <<- 'EOF' > $path/$name_final.$nchr.chr.snape.out.R
x <- read.table(file="'$path/$name_final.$nchr.chr.snape.out'", na.strings="NA", fill=TRUE, header=F)
x.filtered  <- x[as.numeric(x[, 9]) >= 0.90,  ]
write.table(x.filtered, file="'$path/$name_final.snape.filtered.out'", sep="\t", quote=F, col.names=F, row.names=F)
EOF
R --vanilla < $path/$name_final.$nchr.chr.snape.out.R
rm $path/$name_final.$nchr.chr.snape.out.R

