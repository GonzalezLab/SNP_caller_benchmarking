#!/bin/bash

module load SAMtools/1.6-foss-2016b
module load Python/2.7.12-foss-2016b 

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

POOLSNP_PATH=$DIRDATA/PoolSNP
NINDV=$(echo $PLOIDY/2 | bc)
MINFR=$(echo  "print 0.50/$NINDV" | perl)

#Generate mpileup files
for (( REP=1; REP<=$NUMREPS; REP++ ))
do
  file=$DIRDATA/BAM_files/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL1.replicate$REP.haploidcov$AVCOVREF.bam     
  pool1=$(echo $file)
  pool2=$(echo $file | sed 's/POOL1/POOL2/')
  pool3=$(echo $file | sed 's/POOL1/POOL3/')
  out_p1=$(echo $file | sed 's/POOL1/POOLS/')
  out_p2=$(echo $out_p1 | sed 's/BAM_files/MPILEUP_POOLSNP/')
  out=${out_p2%.*}
  samtools mpileup -B -f  $FASTA/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa  $pool1 $pool2 $pool3 > $out.mpileup
done

#run PoolSNP with miss-fraction = 0.8
for (( REP=1; REP<=$NUMREPS; REP++ ))
do
    $POOLSNP_PATH/PoolSNP.sh \
	mpileup=$MPILEUP/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOLS.replicate$REP.haploidcov$AVCOVREF.mpileup \
	output=$POOLSNP/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOLS.replicate$REP.haploidcov$AVCOVREF.mf.output \
	reference=$FASTA/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa \
	names=africa,america,europe \
	min-cov=8 \
	max-cov=0.95 \
	min-count=2 \
	min-freq=$MINFR \
	miss-frac=0.8 \
	bq=15 \
    jobs=1 \
	BS=1
done
