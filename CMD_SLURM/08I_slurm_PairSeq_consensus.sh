#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=ilePairCons # Job name
#SBATCH --partition=hugemem # Partition
#SBATCH --cpus-per-task=50
#SBATCH --output=ILEUM_PairConsensus.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

parallel -j 50 :::: ../CMD_FOF/cmd_pair_Consensus.fof
