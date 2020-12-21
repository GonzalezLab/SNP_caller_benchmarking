#!/bin/bash

module load VarScan/2.4.2-Java-1.8.0_92 
module load Perl/6-Rakudo-star-2018.01
module load R/3.4.2-foss-2016b

DIRDATA=/path_to/Callers
RUNFOLDER=/path_to/Callers/run
FASTQ=/path_to/Callers/FASTQ_files
BAM=/path_to/Callers/BAM_files
PILEUP=/path_to/Callers/PILEUP_files
MPILEUP=/path_to/Callers/MPILEUP_files
FASTA=/path_to/Callers/FASTA_files
MS=/path_to/Callers/MS_files
TRUE=/path_to/Callers/TRUEDATA_results
STATS=/path_to/Callers/STATS_data
MAPGD=/path_to/Callers/MAPGD_results
SNAPE=/path_to/Callers/SNAPE_results
VARSCAN2=/path_to/Callers/VARSCAN2_results
POOLSNP=/path_to/Callers/POOLSNP_results
RUNPREFIX=ms_sim

PLOIDY=$1  #Number of sequences per pool (50 diploid indv -> 100 sequences)
AVCOV=$2   #Coverage per pool 
NUMINV=$3  #Number of individuals per pool 
NUMSEQS=$4 #Total number of sequences summing all pools (i.e, if three pools of 100 sequences each pool = 300)
NUMREPS=$5 #Number of ms replicates to be simulated
THETA=$6   #theta per region (4Nmu*SEQLENGTH)
RHO=$7     #Rho per region (4Nmu*r)
SEQLEN=$8  #Length of the region to be simulated
AVCOVREF=`echo "scale=2; $AVCOV/$PLOIDY" | bc` #Coverage per chromosome
POOLNUM=$(($NUMSEQS/$PLOIDY)) #Number of different pools 

NINDV=$(echo $PLOIDY/2 | bc)
MINFR=$(echo  "print 0.50/$NINDV" | perl)
MINC=8
MINRD=2
PVAL=0.05

#Run VarScan
for (( REP=1; REP<=$NUMREPS; REP++ ))
do
	#file=$DIRDATA/BAM_files/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL1.replicate$REP.haploidcov$AVCOVREF.sorted.bam    
    #pool1=$(echo $file)
    #pool2=$(echo $file | sed 's/POOL1/POOL2/')
    #pool3=$(echo $file | sed 's/POOL1/POOL3/')
    #out_p1=$(echo $file | sed 's/POOL1/POOLS/')
    #out_p2=$(echo $out_p1 | sed 's/BAM_files/MPILEUP_files/')
    #out=${out_p2%.sorted.*}

    #samtools mpileup -f -B $FASTA/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa $pool1 $pool2 $pool3 > $out.mpileup	

	java -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar mpileup2snp $MPILEUP/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOLS.replicate$REP.haploidcov$AVCOVREF.mpileup \
		--min-coverage $MINC --min-reads2 $MINRD --min-var-freq $MINFR --p-value $PVAL \
		2> $VARSCAN2/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOLS.replicate$REP.haploidcov$AVCOVREF.sc.txt \
		1> $VARSCAN2/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOLS.replicate$REP.haploidcov$AVCOVREF.varscan.out.txt
done

#Filter those SNPs in VarScan2 output with the cross-sample P-values for calling variants > 0.00001, 0.0001 and 0.01 and 0.05
for (( REP=1; REP<=$NUMREPS; REP++ ))
do
	cat <<- 'EOF' > var.R    
	args <- commandArgs()
	print(args)
	options(digits=9)
	file <- read.table(file=args[9], header=T, sep="\t")
	subset <- gen <- pop1 <- pop2 <- pop3 <- all <- NULL
	for (line in 1:length(file[, 1])) {
		subset <- rbind(subset, file[line, c(2:4)])
		gen <- rbind(gen, as.matrix(t(unlist(strsplit(as.character(file[line, 5]),":")))))
		pop1 <- rbind(pop1, as.matrix(t(unlist(strsplit(unlist(strsplit(as.character(file[line, 11])," "))[1], ":")))))
		pop2 <- rbind(pop2, as.matrix(t(unlist(strsplit(unlist(strsplit(as.character(file[line, 11])," "))[2], ":")))))
		pop3 <- rbind(pop3, as.matrix(t(unlist(strsplit(unlist(strsplit(as.character(file[line, 11])," "))[3], ":"))))) }
		all <- cbind(subset, gen, pop1, pop2, pop3)
		
		colnames(all) <- c("Position", "Ref", "Alt", "Cons_all", "Cov_all", "Reads_Ref_all", "Reads_Alt_all", "Freq_all", "Pval_all", "Cons_pop1", "Cov_pop1", "Reads_Ref_pop1", "Reads_Alt_pop1", "Freq_pop1", "Pval_pop1",  "Cons_pop2", "Cov_pop2", "Reads_Ref_pop2", "Reads_Alt_pop2", "Freq_pop2", "Pval_pop2", "Cons_pop3", "Cov_pop3", "Reads_Ref_pop3", "Reads_Alt_pop3", "Freq_pop3", "Pval_pop3")
		all$Cov_all  <- as.numeric(as.character(all$Cov_all))
		all$Reads_Ref_all  <- as.numeric(as.character(all$Reads_Ref_all))
		all$Reads_Alt_all  <- as.numeric(as.character(all$Reads_Alt_all))
		all$Pval_all  <- as.numeric(as.character(all$Pval_all))
		all$Cov_pop1  <- as.numeric(as.character(all$Cov_pop1))
		all$Reads_Ref_pop1  <- as.numeric(as.character(all$Reads_Ref_pop1))
		all$Reads_Alt_pop1  <- as.numeric(as.character(all$Reads_Alt_pop1))
		all$Pval_pop1  <- as.numeric(as.character(all$Pval_pop1))
		all$Cov_pop2  <- as.numeric(as.character(all$Cov_pop2))
		all$Reads_Ref_pop2  <- as.numeric(as.character(all$Reads_Ref_pop2))
		all$Reads_Alt_pop2  <- as.numeric(as.character(all$Reads_Alt_pop2))
		all$Pval_pop2  <- as.numeric(as.character(all$Pval_pop2))
		all$Cov_pop3  <- as.numeric(as.character(all$Cov_pop3))
		all$Reads_Ref_pop3  <- as.numeric(as.character(all$Reads_Ref_pop3))
		all$Reads_Alt_pop3  <- as.numeric(as.character(all$Reads_Alt_pop3))
		all$Pval_pop3  <- as.numeric(as.character(all$Pval_pop3))

		all[ ,8]  <- as.matrix(sapply(all[ ,8], gsub,pattern="%",replacement=""))
		all[ ,8] <- as.numeric(all[ ,8])/100
		all[ ,14]  <- as.matrix(sapply(all[ ,14], gsub,pattern="%",replacement=""))
		all[ ,14] <- as.numeric(all[ ,14])/100
		all[ ,20]  <- as.matrix(sapply(all[ ,20], gsub,pattern="%",replacement=""))
		all[ ,20] <- as.numeric(all[ ,20])/100
		all[ ,26]  <- as.matrix(sapply(all[ ,26], gsub,pattern="%",replacement=""))
		all[ ,26] <- as.numeric(all[ ,26])/100

		alpha_22 <- all[which(all[,9] < as.numeric(args[6])), ]
		alpha_11 <- all[which(all[,9] < as.numeric(args[7])), ]
		alpha_6 <- all[which(all[,9] < as.numeric(args[8])), ]

		write.table(alpha_22, file=args[10], sep="\t", quote=F, row.names=F, col.names=T, append=F)
		write.table(alpha_11, file=args[11], sep="\t", quote=F, row.names=F, col.names=T, append=F)
		write.table(alpha_6, file=args[12], sep="\t", quote=F, row.names=F, col.names=T, append=F)
		write.table(all, file=args[13], sep="\t", quote=F, row.names=F, col.names=T, append=F)
		EOF
		Rscript var.R  0.00001 \
		0.001 \
		0.01 \
		$VARSCAN2/$NUMINV.$PLOIDY.$AVCOV/MPILEUP/$RUNPREFIX.POOLS.replicate$REP.haploidcov$AVCOVREF.varscan.out.txt \
		$VARSCAN2/$NUMINV.$PLOIDY.$AVCOV/MPILEUP/$RUNPREFIX.POOLS.replicate$REP.haploidcov$AVCOVREF.varscan.Pval.0.00001.txt \
		$VARSCAN2/$NUMINV.$PLOIDY.$AVCOV/MPILEUP/$RUNPREFIX.POOLS.replicate$REP.haploidcov$AVCOVREF.varscan.Pval.0.001.txt \
		$VARSCAN2/$NUMINV.$PLOIDY.$AVCOV/MPILEUP/$RUNPREFIX.POOLS.replicate$REP.haploidcov$AVCOVREF.varscan.Pval.0.01.txt \
		$VARSCAN2/$NUMINV.$PLOIDY.$AVCOV/MPILEUP/$RUNPREFIX.POOLS.replicate$REP.haploidcov$AVCOVREF.varscan.No_lim_Pval.txt 
		rm var.R
done

