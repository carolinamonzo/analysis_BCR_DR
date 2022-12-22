#!/bin/bash

#SBATCH --mem-per-cpu=10gb #Reserving memory
#SBATCH --job-name=SPL1umi # Job name
#SBATCH --partition=long # Partition
#SBATCH --cpus-per-task=60 # Run on one node
#SBATCH --nodes=2
#SBATCH --output=SPLEEN_R1_Mask_primers_UMI.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

parallel -j 2 :::: ../CMD_FOF/cmd_SplR1GetUMI.fof
