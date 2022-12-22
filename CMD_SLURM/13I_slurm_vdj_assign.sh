#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=VDJile # Job name
#SBATCH --partition=hugemem # Partition
#SBATCH --cpus-per-task=30 # Run on one node
#SBATCH --output=ILEUM_vdj_assign.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

module load blast+
source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

parallel -j 1 :::: ../CMD_FOF/cmd_vdj_assign.fof
