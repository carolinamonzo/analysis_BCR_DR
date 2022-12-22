#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=SQLspl # Job name
#SBATCH --partition=himem # Partition
#SBATCH --nodes=1 # Run on one node
#SBATCH --output=./logs/SPL_productive.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc
module load blast+

parallel -j 20 :::: ../CMD_FOF/cmd_SPL_productive.fof
