#!/bin/bash

#SBATCH --mem-per-cpu=40gb #Reserving memory
#SBATCH --job-name=SPL2split # Job name
#SBATCH --partition=long # Partition
#SBATCH --nodes=5 # Run on one node
#SBATCH --output=SPLEEN_R2_SplitSeq_parallel_slurm.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

parallel -j 5 :::: ../CMD_FOF/cmd_SPLR2SplitSeq.fof
