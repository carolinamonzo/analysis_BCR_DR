#!/bin/bash

#SBATCH --mem-per-cpu=10gb #Reserving memory
#SBATCH --job-name=ILEParFilterSeq # Job name
#SBATCH --partition=blade # Partition
#SBATCH --nodes=5 # Run on one node
#SBATCH --output=ILEUM_Mask_primers_parallel_slurm.log # Standard output and error log

source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

parallel -j 5 :::: ../CMD_FOF/cmd_ILEUMMaskPrimers_BC.fof
