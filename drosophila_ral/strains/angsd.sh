#!/bin/bash

module load angsd/0.925-foss-2016b 

path= #full path to bam files
ref= #full path to reference genome

#filelist is the file with all the name of the bam file files for each of the 30 strains 

#run angsd with pval 1e-4 
angsd -bam $path/bam.filelist -nInd 30 -GL 2 -doGeno 2 -ref $ref -doMajorMinor 4 -doMaf 1 -C 50 -baq 1 -doPost 1 -indF testF_seed1_pval1.indF -SNP_pval 1e-4  -doVcf 1 -out angsd_out_ral_pval1

#run angsd with pval 1e-6 
angsd -bam $path/bam.filelist -nInd 30 -GL 2 -doGeno 2 -ref $ref -doMajorMinor 4 -doMaf 1 -C 50 -baq 1 -doPost 1 -indF testF_seed1_pval2.indF -SNP_pval 1e-6  -doVcf 1 -out angsd_out_ral_pval2
