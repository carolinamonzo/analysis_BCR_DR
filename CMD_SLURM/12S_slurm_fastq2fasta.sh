#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=splfq2fa # Job name
#SBATCH --partition=blade # Partition
#SBATCH --cpus-per-task=10
#SBATCH --output=./logs/SPL_fastq2fasta.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

module load blast+
source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

parallel -j 20 :::: ../CMD_FOF/cmd_SPL_fastq2fasta.fof
