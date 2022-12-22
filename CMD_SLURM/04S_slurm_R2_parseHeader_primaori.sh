#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=SPLprima # Job name
#SBATCH --partition=blade,himem # Partition
#SBATCH --nodes=10 # Run on one node
#SBATCH --output=SPLEEN_ParseHeader_primaori2.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

parallel -j 10 :::: ../CMD_FOF/cmd_SPL_ParseHeader_primaori2.fof
