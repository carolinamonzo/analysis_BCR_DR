#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=ILE1UMIMaskPrimers # Job name
#SBATCH --partition=blade,himem # Partition
#SBATCH --nodes=10 # Run on one node
#SBATCH --output=ILEUM_R1_UMI_Mask_primers_parallel_slurm-2.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc
module load cdhit

#parallel -j 10 :::: ../CMD_FOF/cmd_IleR1_5primaori1.fof
#parallel -j 10 :::: ../CMD_FOF/cmd_IleR1_3primaori1.fof
parallel -j 10 :::: ../CMD_FOF/cmd_IleR1GetUMI.fof
