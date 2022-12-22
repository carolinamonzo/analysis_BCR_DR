#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=SPLiso # Job name
#SBATCH --partition=blade # Partition
#SBATCH --nodes=10
#SBATCH --output=./logs/SPL_parseTableIsotypes.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

module load blast+
source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

parallel -j 10 :::: ../CMD_FOF/cmd_SPL_parseTableIsotypes.fof
