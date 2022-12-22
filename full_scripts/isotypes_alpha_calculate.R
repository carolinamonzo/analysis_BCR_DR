## Isotypes alpha calculate values

setwd("~/CM_IGseq/scripts/full_scripts")

suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(tigger))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(scoper))
#library(cowplot)

remove(list = ls())

text = "ILEUM"

#db <- read_airr("../../analysis/created_germlines/reannotated_db/merged_changeo_ILE.tsv")

db <- read_airr("../../analysis/created_germlines/reannotated_db/merged_changeo_ILE_isosum.tsv")
list_files <- unique(db[["sample_id"]])

for (i in 1:length(list_files)){
  A_iso <- alphaDiversity(db %>% filter(sample_id == list_files[[i]]), group = "sum_iso", clone = "clone_id", min_q=0, max_q=3, step_q=0.1, ci=0.95, nboot=200)
  write.table(A_iso@diversity, file = paste0("../../analysis/plots/alpha_isotypes/alpha_iso", list_files[[i]], "_ILE_isosum.tsv"), append = TRUE, col.names = TRUE, row.names = TRUE)
  write.table(A_iso@tests, file = paste0("../../analysis/plots/alpha_isotypes/alpha_iso", list_files[[i]], "_stats_ILE_isosum.tsv"), append=TRUE, col.names = TRUE, row.names = TRUE)
}

# Clear environment and repeat with spleen

remove(list = ls())

text = "SPLEEN"

db <- read_airr("../../analysis/created_germlines/reannotated_db/merged_changeo_SPL_isosum.tsv")

list_files <- unique(db[["sample_id"]])

for (i in 1:length(list_files)){
  A_iso <- alphaDiversity(db %>% filter(sample_id == list_files[[i]]), group = "sum_iso", clone = "clone_id", min_q=0, max_q=3, step_q=0.1, ci=0.95, nboot=200)
  write.table(A_iso@diversity, file = paste0("../../analysis/plots/alpha_isotypes/alpha_iso", list_files[[i]], "_SPL_isosum.tsv"), append = TRUE, col.names = TRUE, row.names = TRUE)
  write.table(A_iso@tests, file = paste0("../../analysis/plots/alpha_isotypes/alpha_iso", list_files[[i]], "_stats_SPL_isosum.tsv"), append=TRUE, col.names = TRUE, row.names = TRUE)
}
