#!/bin/bash

#SBATCH --mem-per-cpu=5gb #Reserving memory
#SBATCH --job-name=ileR2Consensus # Job name
#SBATCH --partition=himem # Partition
#SBATCH --cpus-per-task=60 # Run on one node
#SBATCH --output=ILEUMR2_umi_consensus.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19249_I6_S6_001_PRIMER-S10_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19249_I6_S6_001_PRIMER-S1_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19249_I6_S6_001_PRIMER-S2_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19249_I6_S6_001_PRIMER-S3_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19249_I6_S6_001_PRIMER-S4_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19249_I6_S6_001_PRIMER-S5_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19249_I6_S6_001_PRIMER-S6_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19249_I6_S6_001_PRIMER-S7_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19249_I6_S6_001_PRIMER-S8_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19249_I6_S6_001_PRIMER-S9_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19250_I7_S7_001_PRIMER-S10_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19250_I7_S7_001_PRIMER-S1_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19250_I7_S7_001_PRIMER-S2_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19250_I7_S7_001_PRIMER-S3_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19250_I7_S7_001_PRIMER-S4_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19250_I7_S7_001_PRIMER-S5_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19250_I7_S7_001_PRIMER-S6_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19250_I7_S7_001_PRIMER-S7_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19250_I7_S7_001_PRIMER-S8_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19250_I7_S7_001_PRIMER-S9_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19251_I8_S8_001_PRIMER-S10_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19251_I8_S8_001_PRIMER-S1_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19251_I8_S8_001_PRIMER-S2_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19251_I8_S8_001_PRIMER-S3_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19251_I8_S8_001_PRIMER-S4_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19251_I8_S8_001_PRIMER-S5_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19251_I8_S8_001_PRIMER-S6_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19251_I8_S8_001_PRIMER-S7_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19251_I8_S8_001_PRIMER-S8_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19251_I8_S8_001_PRIMER-S9_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19252_I9_S9_001_PRIMER-S10_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19252_I9_S9_001_PRIMER-S1_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19252_I9_S9_001_PRIMER-S2_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19252_I9_S9_001_PRIMER-S3_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19252_I9_S9_001_PRIMER-S4_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19252_I9_S9_001_PRIMER-S5_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19252_I9_S9_001_PRIMER-S6_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19252_I9_S9_001_PRIMER-S7_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19252_I9_S9_001_PRIMER-S8_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19252_I9_S9_001_PRIMER-S9_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19253_I10_S10_001_PRIMER-S10_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19253_I10_S10_001_PRIMER-S1_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19253_I10_S10_001_PRIMER-S2_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19253_I10_S10_001_PRIMER-S3_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19253_I10_S10_001_PRIMER-S4_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19253_I10_S10_001_PRIMER-S5_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19253_I10_S10_001_PRIMER-S6_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19253_I10_S10_001_PRIMER-S7_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19253_I10_S10_001_PRIMER-S8_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
BuildConsensus.py -s ../../data/fastq_UMI/rmdup/dedup-mpg_L19253_I10_S10_001_PRIMER-S9_R2_pair-pass.fastq --bf BARCODE --cf PRIMER --nproc 60 --act set --maxerror 0.1  --outdir ../../data/fastq_UMI/build_consensus/
