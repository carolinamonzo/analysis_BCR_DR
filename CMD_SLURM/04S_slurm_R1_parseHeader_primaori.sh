#!/bin/bash

#SBATCH --mem-per-cpu=10gb #Reserving memory
#SBATCH --job-name=SPL1prima # Job name
#SBATCH --partition=blade,himem,hugemem # Partition
#SBATCH --nodes=20 # Run on one node
#SBATCH --output=SPLEEN_ParseHeader_primaori1.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

source /beegfs/group_lp/home/CMonzo/.bashrc
source ~/myconda.sh
conda activate MPI_cmc

parallel -j 20 :::: ../CMD_FOF/cmd_SPL_ParseHeader_primaori1.fof
