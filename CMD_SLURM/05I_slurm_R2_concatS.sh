#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=concatPrimer # Job name
#SBATCH --partition=himem # Partition
#SBATCH --nodes=30 # Run on one node
#SBATCH --output=ILEUM_concat_primer_S_R.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc
module load cdhit

parallel -j 30 :::: ../CMD_FOF/cmd_concat_PRIMER_S_R2.fof
parallel -j 30 :::: ../CMD_FOF/cmd_concat_PRIMER_S_R1.fof
