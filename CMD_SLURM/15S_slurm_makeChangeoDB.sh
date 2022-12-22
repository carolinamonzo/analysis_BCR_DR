#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=SQLspl # Job name
#SBATCH --partition=blade # Partition
#SBATCH --nodes=10 # Run on one node
#SBATCH --output=./logs/SPL_makeChangeoDB.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc
module load blast+

parallel -j 10 :::: ../CMD_FOF/cmd_SPL_makeChangeoDB.fof
