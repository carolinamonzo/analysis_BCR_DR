#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=ileISO # Job name
#SBATCH --partition=himem,hugemem # Partition
#SBATCH --cpus-per-task=10
#SBATCH --output=ILEUM_parseTableIsotypes.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

module load blast+
source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

parallel -j 10 :::: ../CMD_FOF/cmd_parseTableIsotypes.fof
