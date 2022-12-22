#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=ilePair # Job name
#SBATCH --partition=himem # Partition
#SBATCH --cpus-per-task=40 # Run on one node
#SBATCH --nodes=1
#SBATCH --output=ILEUM_umi_Pairseq.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

parallel -j 40 :::: ../CMD_FOF/cmd_umi_PairSeq.fof 
