#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=ILE2UMIMaskPrimers # Job name
#SBATCH --partition=blade,himem # Partition
#SBATCH --nodes=20 # Run on one node
#SBATCH --output=ILEUM_R2_UMI_Mask_primers_parallel_slurm-2.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc
module load cdhit

#parallel -j 20 :::: ../CMD_FOF/cmd_ParseHeader_5primaori2.fof
#parallel -j 20 :::: ../CMD_FOF/cmd_ParseHeader_3primaori2.fof
parallel -j 20 :::: ../CMD_FOF/cmd_IleR2GetUMI.fof
