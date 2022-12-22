#!/bin/bash

#SBATCH --job-name=ILEclones # Job name
#SBATCH --partition=himem # Partition
#SBATCH --cpus-per-task=60 # Run on one node
#SBATCH --nodes=2
#SBATCH --output=ILE_defineClones.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

parallel -j 2 :::: ../CMD_FOF/cmd_ILE_defineClones.fof
