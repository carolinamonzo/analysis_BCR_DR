#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=VDJspl # Job name
#SBATCH --partition=blade # Partition
#SBATCH --cpus-per-task=30 # Run on one node
#SBATCH --nodes=2
#SBATCH --output=./logs/SPL_vdj_assign.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

module load blast+
source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

parallel -j 2 :::: ../CMD_FOF/cmd_SPL_vdj_assign.fof
