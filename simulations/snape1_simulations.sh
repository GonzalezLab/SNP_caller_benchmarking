#!/bin/bash

module load R/3.4.2-foss-2016b
module load snape-pooled/1 

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

#Remove positions with 0 reads
DIR=$PILEUP/$NUMINV.$PLOIDY.$AVCOV
for file in $DIR/*.pileup
do
    if [ $(cut -f4 $file | grep -c "^0$") != 0  ]
	then
		if [ ! -d $DIR/PILEUP_MOD ] 
	    then
	    	mkdir $DIR/PILEUP_MOD
	    fi
		out=$(basename $file) 
		cat $file | grep -v $'\t0\t' > $DIR/PILEUP_MOD/$out
	fi
done    

#Run SNAPE with flat prior
for ((POOL=1; POOL<=$POOLNUM; POOL++)) 
do   
    for ((REP=1; REP<=$NUMREPS; REP++))
    do
	if [ `ls -1A $DIR/PILEUP_MOD | wc -l` != 0 ]
	then 
	    diff -q $DIR $DIR/PILEUP_MOD | grep $DIR | grep -E "^Only in.*pileup" | sed -n 's/[^:]*: //p' > $DIR/PILEUP_MOD/list.txt
	    while IFS= read -r line
	    do
			cp $DIR/$line  $DIR/PILEUP_MOD 
		done < "$DIR/PILEUP_MOD/list.txt"
	fi  
	rm $DIR/PILEUP_MOD/list.txt
	
	if [ ! -d $DIR/PILEUP_MOD ]
	then
	    snape-pooled -nchr $PLOIDY -theta 0.01 -D 0.10 -fold unfolded  -priortype flat < \
		$DIR/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.pileup \
		> $SNAPE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.snape.out
	else
	    snape-pooled -nchr $PLOIDY -theta 0.01 -D 0.10 -fold unfolded  -priortype flat < \
		$DIR/PILEUP_MOD/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.pileup \
		> $SNAPE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.snape.out
	fi
	done
done

#Filtered out SNPs with a posterior probability < 0.90
for (( REP=1; REP<=$NUMREPS; REP++ ))
do
	for (( POOL=1; POOL<=$POOLNUM; POOL++ ))
	do	
		cat <<- 'EOF' > $SNAPE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.snape.out.R  
		x <- read.table(file="'$SNAPE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.snape.out'", na.strings="NA", fill = TRUE, header=F)
		x.filtered  <- x[as.numeric(x[, 9]) >= 0.90,  ]
		write.table(x.filtered, file="'$SNAPE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.snape.filtered.out'", sep="\t", quote=F, col.names=F, row.names=F)
		EOF
		R --vanilla < $SNAPE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.snape.out.R   
	done
done
rm $SNAPE/$NUMINV.$PLOIDY.$AVCOV/*snape.out.R
