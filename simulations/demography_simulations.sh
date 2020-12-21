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
echo "*** COALESCENT SIMULATIONS UNDER THE DEMOGRAPHIC MODEL ***"
echo "**********************************************************"
echo " "
echo " "

for ((REP=1; REP<=$NUMREPS; REP++))
do
	##coalescent simulations 
    Pipeliner writeRandomSeq -outfile $RUNFOLDER/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa -len $SEQLEN #Generates a sequence with the same length as the ms simulations
    Pipeliner seedMS -folder $RUNFOLDER/$NUMINV.$PLOIDY.$AVCOV/ 

    ##ms command line
    ms $NUMSEQS 1 -t $THETA -r $RHO $SEQLEN -I 3 $PLOIDY $PLOIDY $PLOIDY  -n 2 3.21 -n 3 0.63 -g 2 120703.6 -g 3 545.88 -eg 0.0000726 2 0.0 \
      -es 0.0000726 2 0.15 -ej 0.0000726 4 3 -ej 0.0000726 2 1 -ej 0.00955 3 1 -en 0.01191 1 0.0001246 -en 0.1192 1 8425.97 > $RUNFOLDER/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.ms

    Pipeliner ms2fas -in_ms $RUNFOLDER/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.ms \
      -in_anc $RUNFOLDER/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa \
      -out_fas $RUNFOLDER/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.fas -force refseq

    ##separate individuals from each population 
    POOL=1
    START=1
    END=$(($PLOIDY*2))
    SUM=$(($PLOIDY*2))
    LIN=$(($NUMSEQS*2))
    while [ $END -le $LIN ]
    do
    	sed -n $START,$END"p" $RUNFOLDER/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.fas > $RUNFOLDER/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.fas
      	POOL=$(($POOL+1))
      	START=$(($START+$SUM))
      	END=$(($END+$SUM))
    done

	##individuals sequencing 
    for ((POOL=1; POOL<=$POOLNUM; POOL++))
    do
       	art_illumina -i $RUNFOLDER/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.fas -ss HS20 -na -q -p -m 350 -s 25 -l 100 -f $AVCOVREF -o $RUNFOLDER/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.fa.
    done
    
    mv $RUNFOLDER/$NUMINV.$PLOIDY.$AVCOV/*.fq $FASTQ/$NUMINV.$PLOIDY.$AVCOV
    mv $RUNFOLDER/$NUMINV.$PLOIDY.$AVCOV/*.fas $FASTA/$NUMINV.$PLOIDY.$AVCOV
    mv $RUNFOLDER/$NUMINV.$PLOIDY.$AVCOV/*.fa $FASTA/$NUMINV.$PLOIDY.$AVCOV
    mv $RUNFOLDER/$NUMINV.$PLOIDY.$AVCOV/*.ms $MS/$NUMINV.$PLOIDY.$AVCOV
    rm $RUNFOLDER/$NUMINV.$PLOIDY.$AVCOV/seedms

    ##mapping
    MAPQ=20
    bwa index -a is $FASTA/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa
    samtools faidx $FASTA/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa

	##for each pool
	for ((POOL=1; POOL<=$POOLNUM; POOL++))
    do
    	bwa mem $FASTA/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa \
			$FASTQ/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.fa.1.fq \
     	  	$FASTQ/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.fa.2.fq \
     	     > $BAM/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.unsorted.sam

        samtools view -q $MAPQ -Sb -h $BAM/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.unsorted.sam \
     	  > $BAM/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.unsorted.bam

        samtools sort $BAM/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.unsorted.bam > \
     	  $BAM/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.sorted.bam

        java -Xmx20G -XX:+UseSerialGC -jar $PICARDPATH/picard.jar AddOrReplaceReadGroups \
     	  I=/$DIRDATA/BAM_files/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.sorted.bam \
     	  O=/$DIRDATA/BAM_files/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.bam \
     	  ID=1 LB=1 PL=illumina PU=1 SM="POOL$POOL"

    	##some stats
    	samtools flagstat $BAM/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.bam \
     	  > $BAM/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.mapping.stats
				    
        bedtools genomecov -d -ibam $BAM/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.bam \
     	  -g $FASTA/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa \
     	  > $BAM/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.bed

    	##BAM to pileup 
    	samtools mpileup -B -f $FASTA/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa \
     	  $BAM/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.bam \
     	  > $PILEUP/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.haploidcov$AVCOVREF.pileup
    done
done
