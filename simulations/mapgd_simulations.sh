#!/bin/bash

module load SAMtools/1.6-foss-2016b
module load MAPGD/0.4.26

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

#Built a mpileup for each replicate and run MAPGD with three different thresholds for the log-likelihood ratio (LLR) test (-a flag)
alpha=( 22 11 6 )
for (( REP=1; REP<=$NUMREPS; REP++ ))
do
	file=$DIRDATA/BAM_files/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL1.replicate$REP.haploidcov$AVCOVREF.sorted.bam    
    pool1=$(echo $file)
    pool2=$(echo $file | sed 's/POOL1/POOL2/')
    pool3=$(echo $file | sed 's/POOL1/POOL3/')
    out_p1=$(echo $file | sed 's/POOL1/POOLS/')
    out_p2=$(echo $out_p1 | sed 's/BAM_files/MPILEUP_files/')
    out=${out_p2%.sorted.*}

    samtools mpileup -f -B $FASTA/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa $pool1 $pool2 $pool3 > $out.mpileup	
    samtools view -H $DIRDATA/BAM_files/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL1.replicate$REP.haploidcov$AVCOVREF.sorted.bam \
		> $DIRDATA/MAPGD_results/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL1.replicate$REP.haploidcov$AVCOVREF.header 
  
    mapgd sam2idx -H $DIRDATA/MAPGD_results/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL1.replicate$REP.haploidcov$AVCOVREF.header
	  
    mapgd proview -i $out.mpileup \
		-H $DIRDATA/MAPGD_results/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL1.replicate$REP.haploidcov$AVCOVREF.header \
		-o $DIRDATA/MAPGD_results/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOLS.replicate$REP.haploidcov$AVCOVREF
	
	for i in "${alpha[@]}"
    do
    	if [ ! -d  $DIRDATA/MAPGD_results/$NUMINV.$PLOIDY.$AVCOV/alpha_$i ]
       	then
	   		mkdir $DIRDATA/MAPGD_results/$NUMINV.$PLOIDY.$AVCOV/alpha_$i
       fi
       mapgd pool -i $DIRDATA/MAPGD_results/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOLS.replicate$REP.haploidcov$AVCOVREF.pro \
	   -a $i -o $DIRDATA/MAPGD_results/$NUMINV.$PLOIDY.$AVCOV/alpha_$i/$RUNPREFIX.POOLS.replicate$REP.haploidcov$AVCOVREF.alpha_$i.allelefrequency-filtered
    done
done
