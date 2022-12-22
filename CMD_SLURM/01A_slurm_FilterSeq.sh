#!/bin/bash

#SBATCH --job-name=ParFilterSeq # Job name
#SBATCH --partition=blade # Partition
#SBATCH --nodes=20 # Run on one node
#SBATCH --output=FilterSeq_parallel_slurm.log # Standard output and error log

source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

parallel -j 20 :::: cmd_FilterSeq.cmd
