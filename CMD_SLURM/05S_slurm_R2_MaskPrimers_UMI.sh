#!/bin/bash

#SBATCH --mem-per-cpu=10gb #Reserving memory
#SBATCH --job-name=SPL2umi # Job name
#SBATCH --partition=long # Partition
#SBATCH --cpus-per-task=30 # Run on one node
#SBATCH --output=SPLEEN_R2_Mask_primers_UMI.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

parallel -j 1 :::: ../CMD_FOF/cmd_SplR2GetUMI.fof
