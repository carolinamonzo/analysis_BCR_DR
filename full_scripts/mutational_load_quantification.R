setwd("~/CM_IGseq/scripts/full_scripts")

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(tigger))
suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(scoper))
#library(cowplot)

text = "SPL"

db <- read_airr(paste0("../../analysis/created_germlines/reannotated_db/merged_changeo_", text, ".tsv"))

## Calculate total mutational count, R and S combined, for IMGT_V
db2 <- observedMutations(db, sequenceColumn = "sequence_alignment", germlineColumn = "germline_alignment_d_mask", frequency = T, combine = F, nproc = 32)

subdf <- db2[, c("sequence_id", "sample_id", "biogroup", "mu_freq_seq_r", "mu_freq_seq_s")]
write.table(subdf, file = "../../analysis/results_tables/SHM_frequency_", text, ".tsv", col.names = TRUE, row.names = TRUE)
