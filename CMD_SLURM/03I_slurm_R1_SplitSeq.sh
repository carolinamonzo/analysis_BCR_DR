#!/bin/bash

#SBATCH --mem-per-cpu=20gb #Reserving memory
#SBATCH --job-name=ILE1split # Job name
#SBATCH --partition=long # Partition
#SBATCH --nodes=5 # Run on one node
#SBATCH --output=ILEUM_R1_SplitSeq_parallel_slurm.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cmonzo@age.mpg.de

source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

parallel -j 5 :::: ../CMD_FOF/cmd_ILER1SplitSeq.fof
