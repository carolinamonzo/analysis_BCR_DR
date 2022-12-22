#!/bin/bash

#SBATCH --job-name=SPLclones # Job name
#SBATCH --partition=blade # Partition
#SBATCH --nodes=2
#SBATCH --output=./logs/ALL_merge_changeo.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

source ~/.bashrc
source ~/myconda.sh

parallel -j 2 :::: ../CMD_FOF/cmd_merge_changeo.fof
