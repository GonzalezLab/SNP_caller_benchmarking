#!/bin/bash
                                                                                                                                                                              
module load VarScan/2.4.2-Java-1.8.0_92 
module load Perl/6-Rakudo-star-2018.01
module load R/3.5.1-foss-2018b

path=$1 #path to input
name=$2 #name of pool
nchr=$3 #number of chromosomes in the pool
cov=$4 #coverage

minc=8
minfr=$(echo "print 1/$nchr" | perl)
minrd=2
pval=0.05

java -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar pileup2snp $path/$name_final.mpileup \
	--min-coverage $minc --min-reads2 $minrd --min-var-freq $minfr \
	--p-value $pval 1> $path/$name.$nchr.chr.varscan.out.txt 2> $path/$name.$nchr.chr.sc.txt

#Filter those SNPs in VarScan2 output with the cross-sample P-values for calling variants > 0.00001, 0.0001 and 0.01 and 0.05
cat <<- 'EOF' > $path/$name.var.R  
args <- commandArgs()
print(args)
options(digits=9)
file <- read.table(file=args[9], header=T, sep="\t")
file$VarFreq  <- as.character(file$VarFreq)
file$Pvalue <- as.numeric(file$Pvalue)
file$VarFreq <- as.matrix(sapply(file$VarFreq, gsub,pattern="%",replacement=""))
file$VarFreq <- as.numeric(file$VarFreq)/100
alpha_22 <- file[file[,12] < as.numeric(args[6]), ]
alpha_11 <- file[file[,12] < as.numeric(args[7]), ]
alpha_6 <- file[file[,12] < as.numeric(args[8]), ]
write.table(alpha_22, file=args[10], sep="\t", quote=F, row.names=F, col.names=T, append=F)
write.table(alpha_11, file=args[11], sep="\t", quote=F, row.names=F, col.names=T, append=F)
write.table(alpha_6, file=args[12], sep="\t", quote=F, row.names=F, col.names=T, append=F)
write.table(file, file=args[13], sep="\t", quote=F, row.names=F, col.names=T, append=F)
EOF

Rscript $path/$name.var.R  0.00001 \
0.001 \
0.01 \
$path/$name.60.chr.varscan.out.txt \
$path/$name.60.chr.varscan.Pval.0.00001.out.txt \
$path/$name.60.chr.varscan.Pval.0.001.out.txt \
$path/$name.60.chr.varscan.Pval.0.01.out.txt \
$path/$name.60.chr.varscan.No_lim_Pval.out.txt 
rm -f $path/$name.var.R