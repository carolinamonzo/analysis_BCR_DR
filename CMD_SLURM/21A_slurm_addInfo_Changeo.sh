#!/bin/bash

#SBATCH --job-name=SPLclones # Job name
#SBATCH --partition=blade # Partition
#SBATCH --nodes=10
#SBATCH --output=./logs/ALL_addInfo_changeo.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

source ~/.bashrc
source ~/myconda.sh

parallel -j 10 :::: ../CMD_FOF/cmd_addInfo_changeo.fof
