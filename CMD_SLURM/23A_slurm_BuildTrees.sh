#!/bin/bash

#SBATCH --job-name=Trees # Job name
#SBATCH --partition=hooli # Partition
#SBATCH --nodes=2
#SBATCH --output=./logs/ALL_BuildTrees.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

source ~/.bashrc
source ~/myconda.sh

parallel -j 2 :::: ../CMD_FOF/cmd_BuildTrees.fof
