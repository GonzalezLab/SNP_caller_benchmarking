#!/bin/bash

module load ms/171108-foss-2016b
module load Pipeliner/0.2.0-foss-2016b
module load GSL/2.1-foss-2016b
module load ART/2016.06.05-foss-2016b
module load BWA/0.7.16a-foss-2016b
module load SAMtools/1.6-foss-2016b
module load GATK/3.4-46-Java-1.8.0_92
module load picard/2.8.3
module load BEDTools/2.27.1-foss-2016b
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

echo " "
echo "**********************************************************"
echo "***** COALESCENT SIMULATIONS UNDER THE NEUTRAL MODEL *****"
echo "**********************************************************"
echo " "
echo " "

if [[ ! -d $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV ]]; then mkdir $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV; fi

for ((REP=1; REP<=$NUMREPS; REP++))
do
	##coalescent simulations
    Pipeliner writeRandomSeq -outfile $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa -len $SEQLEN #generates a sequence with the same length as the ms simulations
    Pipeliner seedMS -folder $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/ #just in case you want to repeat the same results

    ##ms command line
    ms $NUMSEQS 1 -t $THETA -r $RHO $SEQLEN -I 3 $PLOIDY $PLOIDY $PLOIDY -ej 4 3 2 -ej 4 2 1 > $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.ms
    
    Pipeliner ms2fas -in_ms $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.ms -in_anc $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa -out_fas $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.fas -force refseq
    
    ##separate individuals from each population
    POOL=1
    START=1
    END=$(($PLOIDY*2))
    SUM=$(($PLOIDY*2))
	LIN=$(($NUMSEQS*2))
    while [ $END -le $LIN ]
	do
      sed -n $START,$END"p" $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.fas > $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.fas
      POOL=$(($POOL+1))
      START=$(($START+$SUM))
      END=$(($END+$SUM))
    done
    
    ##individuals sequencing  the variable $POOL is here the individual with ploidy P, in my case is the population with pool of P chromosomes
    for ((POOL=1; POOL<=$POOLNUM; POOL++))
    do
    	art_illumina -i $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.fas -ss HS20 -na -q -p -m 350 -s 25 -l 100 -f $AVCOVREF -o $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.fa.
	done
    
    mv $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/*.fq $FASTQ/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV
    mv $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/*.fas $FASTA/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV
    mv $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/*.fa $FASTA/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV
    mv $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/*.ms $MS/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV
    rm $RUNFOLDER/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/seedms
    
    ##mapping
    MAPQ=20
    bwa index -a is $FASTA/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa
    samtools faidx $FASTA/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa
    
    ##for each pool
    for ((POOL=1; POOL<=$POOLNUM; POOL++)) 
    do
    	bwa mem $FASTA/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa \
     	    $FASTQ/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.fa.1.fq \
     	    $FASTQ/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.fa.2.fq \
     	    > $BAM/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.unsorted.sam
     	samtools view -q $MAPQ -Sb $BAM/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.unsorted.sam \
     	    > $BAM/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.unsorted.bam
     	samtools sort $BAM/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.unsorted.bam > \
    	    $BAM/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.sorted.bam

     	java -Xmx20G -XX:+UseSerialGC -jar $PICARDPATH/picard.jar AddOrReplaceReadGroups I=/$DIRDATA/BAM_files/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.sorted.bam \
     	    O=/$DIRDATA/BAM_files/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.bam ID=1 LB=1 PL=illumina PU=1 SM="POOL$POOL"
	
     	##some stats
     	samtools flagstat $BAM/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.bam \
     	    > $BAM/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.mapping.stats
     	bedtools genomecov -d -ibam $BAM/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.bam \
     	    -g $FASTA/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa > $BAM/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.bed
	
     	##BAM to pileup
     	samtools mpileup -B -f $FASTA/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa $BAM/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.bam \
     	    > $PILEUP/NEUTRAL.$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.pileup
     done
done