#!/bin/bash

module load mstatspop/0.1beta-foss-2016b

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

##run mstatspop for each simulation and pool
for (( POOL=1; POOL<=$POOLNUM; POOL++ )) 
do   
	for (( REP=1; REP<=$NUMREPS; REP++ ))
    do
		cat $FASTA/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.fas > $TRUE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.fas
		cat $FASTA/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.replicate$REP.reference.fa >> $TRUE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.fas
		mstatspop -f fasta -i $TRUE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.fas -u 1 -n refseq -o 10 \
	    -N 2 $PLOIDY 1 -G 1 -T $TRUE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.truedata.results
	done
done

##get SNPs, SFS and SNP frequency from the mstatspop output. 
filecolumns=$TRUE/columns.txt #File with the name of the columns of interest in the mstatspop output
for (( POOL=1; POOL<=$POOLNUM; POOL++ )) 
do   
    for (( REP=1; REP<=$NUMREPS; REP++ ))
    do
		fileinput=$TRUE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.truedata.results
		for line in $(cat $filecolumns) 
		do
	    	if [ $(cat $fileinput | grep -F $line | wc -l) != 0  ]
			then
				if [[ $line =~ "Positions" ]]
		    	then
		    		line2=$(cat $fileinput | grep -n "$line"| awk '{print $1}' FS=":" | tail -1)
		    		line2=$(($line2 + 1))
		    		cat $fileinput | head -$line2 | tail -1 >> $TRUE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.stats.txt
				else                  
		    		echo $(cat $fileinput | grep -F $line) >> $TRUE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.stats.txt
		    	fi
			fi
	    done

		##SNPs
		lin1=$(head -1 $TRUE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.stats.txt)
		sites=( $lin1 )
		for i in "${sites[@]}"
		do
	    	echo $i >> $TRUE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.sites.txt
	    done
	    
		##SFS
		lin2=$(head -2 $TRUE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.stats.txt | tail -1)
		sfs=( $lin2 )
		for ((i = 0; i < ${#sfs[@]}; i+=2))
		do   
	    	echo -e "${sfs[$i]}\t${sfs[$i+1]}" >> $TRUE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.SFS.txt
	    done
	    
		##SNP frequency
		lin3=$(head -3 $TRUE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.stats.txt | tail -1)
		freq=( $lin3 )
		for ((i = 0; i < ${#freq[@]}; i+=2))
		do   
	    	echo -e "${freq[$i]}\t${freq[$i+1]}" >> $TRUE/$NUMINV.$PLOIDY.$AVCOV/$RUNPREFIX.POOL$POOL.replicate$REP.freq.txt
	    done
	done
done