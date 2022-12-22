setwd("~/CM_IGseq/scripts/full_scripts")

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(tigger))
suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(scoper))
#library(cowplot)
# new_day <- gsub("-", "", as.character(Sys.Date()))
# text = "SPL"
# 
# # Read airr
# db <- read_airr(paste0("../../analysis/created_germlines/reannotated_db/merged_changeo_", text,".tsv"))
# 
# ## Now calculating by timepoint
# 
# print("--- Collapsing 5M ---")
# 
# db5 <- subset(db, biogroup %in% c("5DR_lifelong", "5AL_lifelong"))
# 
# clones5 <- collapseClones(db5, cloneColumn = "clone_id", sequenceColumn = "sequence_alignment", germlineColumn = "germline_alignment_d_mask",
#                           method = "thresholdedFreq", minimumFrequency = 0.6, includeAmbiguous = FALSE,
#                           breakTiesStochastic = FALSE, nproc = 32, fields = c("sample_id"))
# 
# saveRDS(clones5, paste0("../../analysis/results_tables/5Mconsensus_clones_selection_", text, "_", new_day, ".rds"))
# 
# print("--- Collapsing 20M ---")
# 
# db20 <- subset(db, biogroup %in% c("20DR_lifelong", "20AL_lifelong", "20AL_DR16M"))
# 
# clones20 <- collapseClones(db20, cloneColumn = "clone_id", sequenceColumn = "sequence_alignment", germlineColumn = "germline_alignment_d_mask",
#                            method = "thresholdedFreq", minimumFrequency = 0.6, includeAmbiguous = FALSE,
#                            breakTiesStochastic = FALSE, nproc = 32, fields = c("sample_id"))
# 
# saveRDS(clones20, paste0("../../analysis/results_tables/20Mconsensus_clones_selection_", text, "_", new_day, ".rds"))
# 
# print("--- Collapsing 24M ---")
# 
# db24 <- subset(db, biogroup %in% c("24DR_lifelong", "24AL_lifelong", "24AL_DR16M", "24AL_DR20M"))
# 
# clones24 <- collapseClones(db24, cloneColumn = "clone_id", sequenceColumn = "sequence_alignment", germlineColumn = "germline_alignment_d_mask",
#                            method = "thresholdedFreq", minimumFrequency = 0.6, includeAmbiguous = FALSE,
#                            breakTiesStochastic = FALSE, nproc = 32, fields = c("sample_id"))
# 
# saveRDS(clones24, paste0("../../analysis/results_tables/24Mconsensus_clones_selection_", text, "_", new_day, ".rds"))
# 
# ### Getting it also by diet through age

print("--- Collapsing AL ---")

remove(list = ls())
text = "SPL"
new_day <- gsub("-", "", as.character(Sys.Date()))
db <- read_airr(paste0("../../analysis/created_germlines/reannotated_db/merged_changeo_", text,".tsv"))

# Constructing clonal consensus sequences, because individual sequences within clonal groups are not really independent events
dbAL <- subset(db, biogroup %in% c("5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong"))

clonesAL <- collapseClones(dbAL, cloneColumn = "clone_id", sequenceColumn = "sequence_alignment", germlineColumn = "germline_alignment_d_mask",
                         method = "thresholdedFreq", minimumFrequency = 0.6, includeAmbiguous = FALSE,
                         breakTiesStochastic = FALSE, nproc = 32, fields = c("sample_id"))

saveRDS(clonesAL, paste0("../../analysis/results_tables/ALconsensus_clones_selection_", text, "_", new_day, ".rds"))

# print("--- Collapsing DR ---")
# 
# # Clean environment and start with next diet
# remove(list = ls())
# text = "SPL"
# new_day <- gsub("-", "", as.character(Sys.Date()))
# db <- read_airr(paste0("../../analysis/created_germlines/reannotated_db/merged_changeo_", text,".tsv"))
# 
# dbDR <- subset(db, biogroup %in% c("5DR_lifelong", "20DR_lifelong", "24DR_lifelong"))
# 
# clonesDR <- collapseClones(dbDR, cloneColumn = "clone_id", sequenceColumn = "sequence_alignment", germlineColumn = "germline_alignment_d_mask",
#                            method = "thresholdedFreq", minimumFrequency = 0.6, includeAmbiguous = FALSE,
#                            breakTiesStochastic = FALSE, nproc = 32, fields = c("sample_id"))
# 
# saveRDS(clonesDR, paste0("../../analysis/results_tables/DRconsensus_clones_selection_", text, "_", new_day, ".rds"))
# 
# print("--- Collapsing AL_DR16M ---")
# # Clean environment and start with next diet
# remove(list = ls())
# text = "SPL"
# new_day <- gsub("-", "", as.character(Sys.Date()))
# db <- read_airr(paste0("../../analysis/created_germlines/reannotated_db/merged_changeo_", text,".tsv"))
# 
# dbALdr16 <- subset(db, biogroup %in% c("5AL_lifelong", "16AL_lifelong", "20AL_DR16M", "24AL_DR16M"))
# 
# clonesALdr16 <- collapseClones(dbALdr16, cloneColumn = "clone_id", sequenceColumn = "sequence_alignment", germlineColumn = "germline_alignment_d_mask",
#                            method = "thresholdedFreq", minimumFrequency = 0.6, includeAmbiguous = FALSE,
#                            breakTiesStochastic = FALSE, nproc = 32, fields = c("sample_id"))
# 
# saveRDS(clonesALdr16, paste0("../../analysis/results_tables/ALDR16Mconsensus_clones_selection_", text, "_", new_day, ".rds"))

print("--- Collapsing AL_DR20M ---")
# Clean environment and start with next diet
remove(list = ls())
text = "SPL"
new_day <- gsub("-", "", as.character(Sys.Date()))
db <- read_airr(paste0("../../analysis/created_germlines/reannotated_db/merged_changeo_", text,".tsv"))

dbALdr20 <- subset(db, biogroup %in% c("5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_DR20M"))

clonesALdr20 <- collapseClones(dbALdr20, cloneColumn = "clone_id", sequenceColumn = "sequence_alignment", germlineColumn = "germline_alignment_d_mask",
                               method = "thresholdedFreq", minimumFrequency = 0.6, includeAmbiguous = FALSE,
                               breakTiesStochastic = FALSE, nproc = 32, fields = "sample_id")

saveRDS(clonesALdr20, paste0("../../analysis/results_tables/ALDR20Mconsensus_clones_selection_", text, "_", new_day, ".rds"))


# ## Getting out stats and data for plotting
# 
# li = c("5M", "20M", "24M", "AL", "DR", "ALDR16M", "ALDR20M")
# 
# for (i in 1:length(li)){
#   
#   gr <- li[[i]]
#   print(gr)
#   
#   
#   clones <- readRDS(paste0("../analysis/results_tables/", gr, "consensus_clones_selection_", text, "_20210927.rds"))
#   
#   # Calculate selection scores from scratch
#   baseline_sub <- calcBaseline(clones, testStatistic="focused",  nproc=32)
#   
#   # Group by sample ID and isotype
#   subject_grouped <- groupBaseline(baseline_sub, groupBy = c("sample_id", "biogroup", "isotype"))
#   
#   # Group the output by biogroup
#   status_grouped <- groupBaseline(subject_grouped, groupBy = c("biogroup"))
# 
#   pdf(paste0("../analysis/plots/selection_density/", gr, "_density_", text, "-", new_day, ".pdf"), width = 6, height = 1.5)
#   plotBaselineDensity(status_grouped, "biogroup", colorValues = sample_colors, sigmaLimits = c(-1,1))
#   dev.off()
#   
#   print(testBaseline(status_grouped, groupBy="biogroup"))
#   
#   write.table(status_grouped@stats, file= paste0("../analysis/results_tables/", gr, "stats_clones_selection_", text, "_20210928.tsv"), col.names = TRUE, row.names = TRUE)
#   
# }