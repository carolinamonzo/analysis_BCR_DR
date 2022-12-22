#!/bin/bash

#SBATCH --mem-per-cpu=10gb #Reserving memory
#SBATCH --job-name=SPL2barcode # Job name
#SBATCH --partition=long # Partition
#SBATCH --cpus-per-task=60 # Run on one node
#SBATCH --output=SPLEEN_R2_Mask_primers_parallel_slurm.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

###parallel -j 5 :::: ../CMD_FOF/cmd_R2_SPLEENMaskPrimers_BC.fof
MaskPrimers.py align -s ../../data/fastq_filtered/mpg_L19244_S1_S1_R2_001_quality-pass.fastq --nproc 60 --outdir ../../data/fastq_barcode/annotated_BC/ --bf BARCODE -p ../../metadata/presto_PCR2_BC.fasta --mode cut
MaskPrimers.py align -s ../../data/fastq_filtered/mpg_L19245_S2_S2_R2_001_quality-pass.fastq --nproc 60 --outdir ../../data/fastq_barcode/annotated_BC/ --bf BARCODE -p ../../metadata/presto_PCR2_BC.fasta --mode cut
MaskPrimers.py align -s ../../data/fastq_filtered/mpg_L19246_S3_S3_R2_001_quality-pass.fastq --nproc 60 --outdir ../../data/fastq_barcode/annotated_BC/ --bf BARCODE -p ../../metadata/presto_PCR2_BC.fasta --mode cut
MaskPrimers.py align -s ../../data/fastq_filtered/mpg_L19247_S4_S4_R2_001_quality-pass.fastq --nproc 60 --outdir ../../data/fastq_barcode/annotated_BC/ --bf BARCODE -p ../../metadata/presto_PCR2_BC.fasta --mode cut
MaskPrimers.py align -s ../../data/fastq_filtered/mpg_L19248_S5_S5_R2_001_quality-pass.fastq --nproc 60 --outdir ../../data/fastq_barcode/annotated_BC/ --bf BARCODE -p ../../metadata/presto_PCR2_BC.fasta --mode cut
