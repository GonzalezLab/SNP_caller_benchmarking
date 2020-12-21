#!/bin/bash

module load angsd/0.925-foss-2016b 

path= #full path to bam files
ref= #full path to reference genome

#filelist is the file with all the name of the bam file files for each of the 30 strains 

#run angsd with pval 1e-4 
angsd -bam $path/bam.filelist -nInd 30 -GL 2 -doGlf 3 -ref $ref -doMajorMinor 4 -doMaf 1 -SNP_pval 1e-4 -out testF_pval1.HWE 

#run angsd with pval 1e-6
angsd -bam $path/bam.filelist -nInd 30 -GL 2 -doGlf 3 -ref $ref -doMajorMinor 4 -doMaf 1 -SNP_pval 1e-6 -out testF_pval2.HWE 
 