#!/bin/bash 
                                                                                                                                                                             
module load SAMtools/1.6-foss-2016b
module load VarScan/2.4.2-Java-1.8.0_92 
module load Perl/6-Rakudo-star-2018.01
module load R/3.4.2-foss-2016b

path=$1 #path to input
name1=$2 #name of pool1
name2=$3 #name of pool2
nchr=$4 #number of chromosomes
cov=$5 #coverage

minc=8
minfr=$(echo "print 1/$nchr" | perl)
minrd=2
pval=0.05

#Generate mpileup
samtools mpileup -f -B $path/dmel_6.12.fa \
	$path/$name1_final.bam \
	$path/$name2_final.bam \
    > $path/$name1.$name2.mpileup

#run VarScan2
java -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar mpileup2snp $path/$name1.$name2.mpileup \
	--min-coverage $minc --min-reads2 $minrd --min-var-freq $minfr \
	--p-value $pval 1> $path/$name1.$name2.$nchr.chr.varscan.out.txt 2> $path/$name1.$name2.$nchr.chr.sc.txt
	
#Filter those SNPs in VarScan2 output with the cross-sample P-values for calling variants > 0.00001, 0.0001 and 0.01 and 0.05
cat <<- 'EOF' > $path/$name1.$name2.var.R 
args <- commandArgs()
print(args)
options(digits=9)
input <- read.table(file=args[9], header=T, sep="\t")
p1 <- as.data.frame(unlist(lapply(strsplit(as.character(input[[11]]), " ", fixed = T), "[", 1)))
p2 <- as.data.frame(unlist(lapply(strsplit(as.character(input[[11]]), " ", fixed = T), "[", 2)))
input <- cbind(input[,1:10], p1, p2)
input <- input[!grepl("*/", as.character(input[,5])), ]
input <- input[!grepl("*/", as.character(input[,11])) & !grepl("*/", as.character(input[,12])), ]
file <- cbind(as.character(input[,1]),as.numeric(input[,2]), 
	as.character(unlist(lapply(strsplit(as.character(input[[5]]), ":", fixed = T), "[", 5))), 
	as.numeric(unlist(lapply(strsplit(as.character(input[[5]]), ":", fixed = T), "[", 6))), 
	as.character(unlist(lapply(strsplit(as.character(input[[11]]), ":", fixed = T), "[", 5))), 
	as.numeric(unlist(lapply(strsplit(as.character(input[[11]]), ":", fixed = T), "[", 6))), 
	as.character(unlist(lapply(strsplit(as.character(input[[12]]), ":", fixed = T), "[", 5))), 
	as.numeric(unlist(lapply(strsplit(as.character(input[[12]]), ":", fixed = T), "[", 6))))
file <- as.data.frame(file)
file <- na.omit(file)
file[,3] <- as.matrix(sapply(file[,3], gsub,pattern="%",replacement=""))
file[,5] <- as.matrix(sapply(file[,5], gsub,pattern="%",replacement=""))
file[,7] <- as.matrix(sapply(file[,7], gsub,pattern="%",replacement=""))
file[,3] <- as.numeric(as.character(file[,3]))/100
file[,5] <- as.numeric(as.character(file[,5]))/100
file[,7] <- as.numeric(as.character(file[,7]))/100
file[,1] <- as.character(file[,1])
file[,2] <- as.numeric(as.character(file[,2]))
file[,3] <- as.numeric(as.character(file[,3]))
file[,4] <- as.numeric(as.character(file[,4]))
file[,5] <- as.numeric(as.character(file[,5]))
file[,6] <- as.numeric(as.character(file[,6]))
file[,7] <- as.numeric(as.character(file[,7]))
file[,8] <- as.numeric(as.character(file[,8]))
alpha_22.p <- file[file[,4] < as.numeric(args[6]), ]
alpha_11.p <- file[file[,4] < as.numeric(args[7]), ]
alpha_6.p <- file[file[,4] < as.numeric(args[8]), ]
colnames(alpha_22.p) <- c("Chrom","Position","Freq_p1p2","Pval_p1p2","Freq_p1","Pval_p1","Freq_p2","Pval_p2")
colnames(alpha_11.p) <- c("Chrom","Position","Freq_p1p2","Pval_p1p2","Freq_p1","Pval_p1","Freq_p2","Pval_p2")
colnames(alpha_6.p) <- c("Chrom","Position","Freq_p1p2","Pval_p1p2","Freq_p1","Pval_p1","Freq_p2","Pval_p2")
write.table(alpha_22.p, file=args[10], sep="\t", quote=F, row.names=F, col.names=T)
write.table(alpha_11.p, file=args[11], sep="\t", quote=F, row.names=F, col.names=T)
write.table(alpha_6.p, file=args[12], sep="\t", quote=F, row.names=F, col.names=T)
write.table(file, file=args[13], sep="\t", quote=F, row.names=F, col.names=T)
EOF

Rscript $path/$name1.$name2.var.R  0.00001 0.001 0.01 \
$path/$name1.$name2.$nchr.chr.varscan.out.txt \
$path/$name1.$name2.join.pools.$nchr.chr.varscan.Pval.0.00001.out.txt \
$path/$name1.$name2.join.pools.$nchr.chr.varscan.Pval.0.001.out.txt \
$path/$name1.$name2.join.pools.$nchr.chr.varscan.Pval.0.01.out.txt \
$path/$name1.$name2.join.pools.$nchr.chr.varscan.No_lim_Pval.out.txt 
rm -f $path/$name1.$name2.var.R
