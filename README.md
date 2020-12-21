# SNP_caller_benchmarking
simulations folder:
1.demography_simulations.sh	     #For running the simulations under the joint demographic model and getting the bam and the mpileup files
2.neutral_simulations.sh	     #For running the simulations under the SNM and getting the bam and the mpileup files
3.mstatspop_simulations.sh 	     #For computing the Number of SNPs, SFS and the SNP Frequency of simulated data
4.mapgd_simulations.sh		     #For calling SNPs using MAPGD with three different thresholds for the log-likelihood ratio (LLR) test
5.varscan2_simulations.sh	     #For calling SNPs using VarScan2 and filtering the SNPs with different cross-sample P-values
6.snape1_simulations.sh		     #For calling SNPs using SNAPE with flat prior and filtering the SNPs with P > 0.90
7.snape2_simulations.sh	     	     #For calling SNPs using SNAPE with informative prior and filtering the SNPs with P > 0.90 
8.poolsnp1_simulations.sh	     #For calling SNPs using PoolSNP with miss-fraction = 0.1 
9.poolsnp2_simulations.sh	     #For calling SNPs using PoolSNP with miss-fraction = 0.8 

drosophila_ral folder:
1.strains folder:
1.1.mapping_strains.sh		     #For trimming, mapping and getting the bam and the mpileup files for each DGRP strain
1.2.angsd4ngsd.sh  	       	     #For estimate the genotype likelihoods for each site of each DGRP strain to calculate their inbreeding coefficient
1.3.ngsd.sh	     	      	     #For estimate the inbreeding coefficient of each DGRP strain using NGSD
1.4.angsd.sh	     	     	     #For calling biallelic SNPs using ANGSD taking into account the inbreeding coefficients of each strain 
2.pools folder:
2.1.mapping_pool_ral.sh		     #For trimming, mapping and getting the bam and the mpileup files for each pool of DGRP strains
2.2.mapgd_pool_ral.sh      	     #For calling SNPs using MAPGD with three different thresholds for the log-likelihood ratio (LLR) test in each pool of DGRP strains
2.3.varscan2_pool_ral.sh 	     #For calling SNPs using VarScan2 and filtering the SNPs with different cross-sample P-values in each pool of DGRP strains
2.4.varscan2_joint_pool_ral.sh       #For calling SNPs using VarScan2 and filtering the SNPs with different cross-sample P-values in pools of DGRP strains joint variant calling option	
2.5.snape1_pool_ral.sh	     	     #For calling SNPs using SNAPE with flat prior and filtering the SNPs with P > 0.90 in each pool of DGRP strains						
2.6.snape2_pool_ral.sh               #For calling SNPs using SNAPE with informative prior and filtering the SNPs with P > 0.90 in each pool of DGRP strains
2.7.poolsnp_pool_ral.sh	     	     #For calling SNPs using PoolSNP in each pool of DGRP strains
2.8.poolsnp_joint_pools_ral.sh       #For calling SNPs using PoolSNP in pools of DGRP strains joint variant calling option with different values of miss-fraction
