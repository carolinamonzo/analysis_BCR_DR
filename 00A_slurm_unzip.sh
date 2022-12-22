#!/bin/bash

#SBATCH --job-name=FilterSeq # Job name
#SBATCH --partition=blade # Partition
#SBATCH --nodes=5 # Run on one node
#SBATCH --cpus-per-task=4 # Number of cpus
#SBATCH --mem-per-cpu=5gb #Job memory request
#SBATCH --time=2-00:00:00 # Time limit days-hrs:min:sec
#SBATCH --output=unzip_slurm.log # Standard output and error log

source /beegfs/group_lp/home/CMonzo/.bashrc
conda activate MPI_cmc

gunzip -c ../data/fastq_raw/mpg_L19244_S1_S1_R1_001.fastq.gz > ../data/fastq_raw/mpg_L19244_S1_S1_R1_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19244_S1_S1_R2_001.fastq.gz > ../data/fastq_raw/mpg_L19244_S1_S1_R2_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19245_S2_S2_R1_001.fastq.gz > ../data/fastq_raw/mpg_L19245_S2_S2_R1_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19245_S2_S2_R2_001.fastq.gz > ../data/fastq_raw/mpg_L19245_S2_S2_R2_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19246_S3_S3_R1_001.fastq.gz > ../data/fastq_raw/mpg_L19246_S3_S3_R1_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19246_S3_S3_R2_001.fastq.gz > ../data/fastq_raw/mpg_L19246_S3_S3_R2_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19247_S4_S4_R1_001.fastq.gz > ../data/fastq_raw/mpg_L19247_S4_S4_R1_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19247_S4_S4_R2_001.fastq.gz > ../data/fastq_raw/mpg_L19247_S4_S4_R2_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19248_S5_S5_R1_001.fastq.gz > ../data/fastq_raw/mpg_L19248_S5_S5_R1_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19248_S5_S5_R2_001.fastq.gz > ../data/fastq_raw/mpg_L19248_S5_S5_R2_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19249_I6_S6_R1_001.fastq.gz > ../data/fastq_raw/mpg_L19249_I6_S6_R1_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19249_I6_S6_R2_001.fastq.gz > ../data/fastq_raw/mpg_L19249_I6_S6_R2_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19250_I7_S7_R1_001.fastq.gz > ../data/fastq_raw/mpg_L19250_I7_S7_R1_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19250_I7_S7_R2_001.fastq.gz > ../data/fastq_raw/mpg_L19250_I7_S7_R2_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19251_I8_S8_R1_001.fastq.gz > ../data/fastq_raw/mpg_L19251_I8_S8_R1_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19251_I8_S8_R2_001.fastq.gz > ../data/fastq_raw/mpg_L19251_I8_S8_R2_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19252_I9_S9_R1_001.fastq.gz > ../data/fastq_raw/mpg_L19252_I9_S9_R1_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19252_I9_S9_R2_001.fastq.gz > ../data/fastq_raw/mpg_L19252_I9_S9_R2_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19253_I10_S10_R1_001.fastq.gz > ../data/fastq_raw/mpg_L19253_I10_S10_R1_001.fastq
gunzip -c ../data/fastq_raw/mpg_L19253_I10_S10_R2_001.fastq.gz > ../data/fastq_raw/mpg_L19253_I10_S10_R2_001.fastq
