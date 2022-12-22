#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=R2concatPrimer # Job name
#SBATCH --partition=blade # Partition
#SBATCH --nodes=20 # Run on one node
#SBATCH --output=SPLEEN_concat_primer_S_R2.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

parallel -j 20 :::: ../CMD_FOF/cmd_SPL_concat_PRIMER_S_R2.fof
