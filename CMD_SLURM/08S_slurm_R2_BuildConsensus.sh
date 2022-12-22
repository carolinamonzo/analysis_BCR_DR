#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=SPL2Consensus # Job name
#SBATCH --partition=hugemem # Partition
#SBATCH --cpus-per-task=60 # Run on one node
#SBATCH --nodes=1
#SBATCH --output=SPL_umi_ConsensusR2.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

parallel -j 1 :::: ../CMD_FOF/cmd_SPL_umi_BuildConsensusR2.fof
