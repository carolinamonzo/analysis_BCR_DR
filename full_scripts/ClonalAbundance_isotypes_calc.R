setwd("~/CM_IGseq/scripts/full_scripts")

suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(tigger))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(scoper))
library(cowplot)

organ = "ILE"
#organ = "SPL"

db <- read_airr(paste0("../../analysis/created_germlines/reannotated_db/merged_changeo_", organ, ".tsv"))

ab1 <- estimateAbundance(db %>% filter(isotype == "IGA"), group = "sample_id", clone = "clone_id", ci = 0.95, nboot = 200, min_n = 5)
write.table(ab1@abundance, file = paste0("../../analysis/plots/clonal_abundance/clab_IGA_", organ, ".tsv"), append = TRUE, col.names = TRUE, row.names = TRUE)
ab2 <- estimateAbundance(db %>% filter(isotype == "IGM"), group = "sample_id", clone = "clone_id", ci = 0.95, nboot = 200, min_n = 5)
write.table(ab2@abundance, file = paste0("../../analysis/plots/clonal_abundance/clab_IGM_", organ, ".tsv"), append = TRUE, col.names = TRUE, row.names = TRUE)
ab3 <- estimateAbundance(db %>% filter(isotype == "IGD"), group = "sample_id", clone = "clone_id", ci = 0.95, nboot = 200, min_n = 5)
write.table(ab3@abundance, file = paste0("../../analysis/plots/clonal_abundance/clab_IGD_", organ, ".tsv"), append = TRUE, col.names = TRUE, row.names = TRUE)
ab4 <- estimateAbundance(db %>% filter(isotype == "IGE"), group = "sample_id", clone = "clone_id", ci = 0.95, nboot = 200, min_n = 5)
write.table(ab4@abundance, file = paste0("../../analysis/plots/clonal_abundance/clab_IGE_", organ, ".tsv"), append = TRUE, col.names = TRUE, row.names = TRUE)
ab5 <- estimateAbundance(db %>% filter(isotype == "IGG12"), group = "sample_id", clone = "clone_id", ci = 0.95, nboot = 200, min_n = 5)
write.table(ab5@abundance, file = paste0("../../analysis/plots/clonal_abundance/clab_IGG12_", organ, ".tsv"), append = TRUE, col.names = TRUE, row.names = TRUE)
ab6 <- estimateAbundance(db %>% filter(isotype == "IGG3"), group = "sample_id", clone = "clone_id", ci = 0.95, nboot = 200, min_n = 5)
write.table(ab6@abundance, file = paste0("../../analysis/plots/clonal_abundance/clab_IGG3_", organ, ".tsv"), append = TRUE, col.names = TRUE, row.names = TRUE)
