#!/bin/bash

#SBATCH --mem-per-cpu=10gb #Reserving memory
#SBATCH --job-name=SPLpair # Job name
#SBATCH --partition=blade # Partition
#SBATCH --nodes=20 # Run on one node
#SBATCH --output=SPLEEN_pairConcat.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

module load blast+
source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

parallel -j 20 :::: ../CMD_FOF/cmd_SPL_umi_PairSeq.fof
