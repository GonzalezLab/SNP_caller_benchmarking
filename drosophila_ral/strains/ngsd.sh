#!/bin/bash

module load ngsF/1.2.0-foss-2016b

path= #path to input file

#run ngsf with seed 1 and for input testF_pval1.HWE.mafs.gz
N_SITES=$((`zcat testF_pval1.HWE.mafs.gz | wc -l`-1)) 
zcat $path/testF_pval1.HWE.glf.gz | ngsF --n_ind 30 --n_sites $N_SITES --glf - --min_epsilon 1e-9 --out testF_seed1.pval1.indF --seed 1234

#run ngsf with seed 2 and for input testF_pval1.HWE.mafs.gz
N_SITES=$((`zcat testF_pval1.HWE.mafs.gzfs.gz | wc -l`-1))
zcat $path/testF_pval1.HWE.glf.gz | ngsF --n_ind 30 --n_sites $N_SITES --glf - --min_epsilon 1e-9 --out testF_seed2.pval1.indF --seed 1456

#run ngsf with seed 3 and for input testF_pval1.HWE.mafs.gz
N_SITES=$((`zcat testF_pval1.HWE.mafs.gz | wc -l`-1))
zcat $path/testF_pval1.HWE.glf.gz | ngsF --n_ind 30 --n_sites $N_SITES --glf - --min_epsilon 1e-9 --out testF_seed3.pval1.indF --seed 4678


#run ngsf with seed 1 and for input testF_pval2.HWE.mafs.gz
N_SITES=$((`zcat testF_pval2.HWE.mafs.gz | wc -l`-1)) 
zcat $path/testF_pval2.HWE.glf.gz | ngsF --n_ind 30 --n_sites $N_SITES --glf - --min_epsilon 1e-9 --out testF_seed1.pval2.indF --seed 1234

#run ngsf with seed 2 and for input testF_pval2.HWE.mafs.gz
N_SITES=$((`zcat testF_pval2.HWE.mafs.gz | wc -l`-1))
zcat $path/testF_pval2.HWE.glf.gz | ngsF --n_ind 30 --n_sites $N_SITES --glf - --min_epsilon 1e-9 --out testF_seed2.pval2.indF --seed 1456

#run ngsf with seed 3 and for input testF_pval2.HWE.mafs.gz
N_SITES=$((`zcat testF_pval2.HWE.mafs.gz | wc -l`-1))
zcat $path/testF_pval2.HWE.glf.gz | ngsF --n_ind 30 --n_sites $N_SITES --glf - --min_epsilon 1e-9 --out testF_seed3.pval2.indF --seed 4678
