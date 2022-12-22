#!/bin/bash

#SBATCH --job-name=SPLpair # Job name
#SBATCH --partition=blade # Partition
#SBATCH --nodes=20
#SBATCH --output=./logs/SPL_pairConsensus.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

parallel -j 20 :::: ../CMD_FOF/cmd_SPL_pairConsensus.fof
