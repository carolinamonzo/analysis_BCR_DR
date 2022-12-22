suppressPackageStartupMessages(library(mixOmics))
suppressPackageStartupMessages(library(dplyr))
theme_set(theme_bw())
setwd("~/workspace/CM_IgSeq_tmp/mixomics/")
path <- "~/workspace/CM_IgSeq_tmp/mixomics/"

metadata <- read.table("~/workspace/16S_final/CM_16S_cross-sectionalF245/metadata/metadata_F245_16S.csv", sep = ";", header = T, row.names = 1, stringsAsFactors = F, comment.char = "")

ile <- read.table("~/workspace/CM_IgSeq_tmp/mixomics/tumor_allBCRmetrics_20220404.csv", sep = ";", header = T, stringsAsFactors = F)
ls_cols <- c("Animal_ID","q", "d","p20","mu_freq_seq_r","mu_freq_seq_s",
             "mut_status",         "cdr3_mu",            "cdr3_sigma",         "IgM-IgD-", 
             "IgM-IgD-SHM-",       "IgM-IgD-SHM--1",     "RDI_uniqueness",     "IGHV1",
             "IGHV10",             "IGHV11",             "IGHV12",             "IGHV13",
             "IGHV14",             "IGHV15",             "IGHV16",             "IGHV2",
             "IGHV3",              "IGHV4",              "IGHV5",              "IGHV6",
             "IGHV7",              "IGHV8",              "IGHV9",              "IGHJ1",
             "IGHJ2",              "IGHJ3",              "IGHJ4",              "relAb_IGA",
             "relAb_IGD",          "relAb_IGE",          "relAb_IGG",          "relAb_IGM",
             "alpha_IGA",          "alpha_IGD",          "alpha_IGE",          "alpha_IGM",
             "alpha_IGG",          "RDI_uniqueness_IGA", "RDI_uniqueness_IGD", "RDI_uniqueness_IGE",
             "RDI_uniqueness_IGM", "RDI_uniqueness_IGG", "p20_IGA",            "p20_IGM",
             "p20_IGD",            "p20_IGE",            "p20_IGG",            "mu_count_seq_r_IGA",
             "mu_count_seq_s_IGA", "mut_status_IGA",     "mu_count_seq_r_IGM", "mu_count_seq_s_IGM",
             "mut_status_IGM",     "mu_count_seq_r_IGG", "mu_count_seq_s_IGG", "mut_status_IGG",
             "mu_count_seq_r_IGD", "mu_count_seq_s_IGD", "mut_status_IGD",     "mu_count_seq_r_IGE",
             "mu_count_seq_s_IGE", "mut_status_IGE",     "cdr3_mu_IGA",        "cdr3_sigma_IGA",
             "cdr3_mu_IGE",        "cdr3_sigma_IGE",     "cdr3_mu_IGD",        "cdr3_sigma_IGD",
             "cdr3_mu_IGG",        "cdr3_sigma_IGG",     "cdr3_mu_IGM",        "cdr3_sigma_IGM")
ile <- ile[, colnames(ile) %in% ls_cols]
ile_new <- ile[ile$q == 0, ]
ile_new$Richness <- ile_new$d
ile_new$Shannon <- ile[ile$q == 1, ]$d
ile_new$Simpson <- ile[ile$q == 2, ]$d
ile_new$q <- NULL
ile_new$d <- NULL
rownames(ile_new) <- ile_new$Animal_ID
ile_new$Animal_ID <- NULL
ile <- t(ile_new)

spl <- read.table("~/workspace/CM_IgSeq_tmp/mixomics/tumor_allBCRmetrics_ILE_20220412.csv", sep = ";", header = T, stringsAsFactors = F)
spl <- spl[, colnames(spl) %in% ls_cols]
spl_new <- spl[spl$q == 0, ]
spl_new$Richness <- spl_new$d
spl_new$Shannon <- spl[spl$q == 1, ]$d
spl_new$Simpson <- spl[spl$q == 2, ]$d
spl_new$q <- NULL
spl_new$d <- NULL
rownames(spl_new) <- spl_new$Animal_ID
spl_new$Animal_ID <- NULL
spl <- t(spl_new)

microbiome <-  read.table("~/workspace/16S_final/CM_16S_cross-sectionalF245/analysis/seqtab_merge3/mergedQC/norm-CLEAN_ASVs_counts_merged_20210615.tsv", sep = "\t", header = T, row.names = 1, stringsAsFactors = F, comment.char = "")
colnames(microbiome) <- gsub("\\.", "-", colnames(microbiome))

# Keep only samples where we also have igseq
microbiome <- microbiome[, colnames(microbiome) %in% colnames(spl)]
microbiome <- microbiome[, colnames(microbiome) %in% colnames(ile)]
ile <- ile[, colnames(ile) %in% colnames(microbiome)]
ile[is.na(ile)] <- 0
spl <- spl[, colnames(spl) %in% colnames(microbiome)]
spl[is.na(spl)] <- 0

ile <- ile[, order(colnames(ile))]
spl <- spl[, order(colnames(spl))]
microbiome <- microbiome[, order(colnames(microbiome))]



metadata <- metadata[metadata$ID %in% colnames(microbiome), ]
##IF all diets, comment
metadata <- metadata[metadata$Treatment %in% c("AL_lifelong", "DR_lifelong"), ]
metadata <- metadata[metadata$Months %in% c(5, 20, 24)]
ile <- ile[, colnames(ile) %in% metadata$ID]
spl <- spl[, colnames(spl) %in% metadata$ID]
microbiome <- microbiome[, colnames(microbiome) %in% metadata$ID]
#


diet <- metadata[, c("Treatment")]
months <- metadata[, c("Months")]

Y <- diet

# Removing clones only found in one sample as they dont add variability
ile <- t(ile)
presabsile <- ile
presabsile[presabsile>0] <-1
presabsile <- presabsile[, which(colSums(presabsile) > 1)]
ile <-  ile[, colnames(ile) %in% colnames(presabsile)]

spl <- t(spl)
presabspl <- spl
presabspl[presabspl>0] <-1
presabspl <- presabspl[, which(colSums(presabspl) > 1)]
spl <-  spl[, colnames(spl) %in% colnames(presabspl)]

spl<-as.data.frame(spl)
ile<-as.data.frame(ile)
microbiome <- as.data.frame(microbiome)

# Now that we have the same order of columns, transpose and put together
data = list(ile = ile, spl = spl,
            microbiome = t(microbiome))
# check dimension
lapply(data, dim)

design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0

design 

sgccda.res = block.splsda(X = data, Y = Y, ncomp = 5, 
                          design = design)

set.seed(123) # for reproducibility, only when the `cpus' argument is not used
# this code takes a couple of min to run
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 4, nrepeat = 10)

#perf.diablo  # lists the different outputs
pdf(paste0(path,"PCcontributions.pdf"))
plot(perf.diablo) 
dev.off()

perf.diablo$choice.ncomp$WeightedVote

# By looking at the PC plots, this combination is most logical
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.ER", "max.dist"]

#Tuning
test.keepX = list(ile = c(5:9, seq(10, 18, 2), seq(20,30,5)), microbiome = c(5:9, seq(10, 18, 2), seq(20,30,5)), spl = c(5:9, seq(10, 18, 2), seq(20,30,5)))

BPPARAM <- BiocParallel::SnowParam(workers = max(parallel::detectCores()-1, 2))

tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, test.keepX = test.keepX, design = design, validation = 'Mfold', folds = 4, nrepeat = 1, dist = "max.dist", BPPARAM = BPPARAM)

list.keepX = tune.TCGA$choice.keepX
list.keepX

# Final model
sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design)
# List the functions related to the object
sgccda.res

## Y has been included automatically in the design, so that the covariance between each block's component and the outcome is maximised

# Extract the selected variables
varmicrobiome <- selectVar(sgccda.res, block = 'microbiome', comp = 1)$microbiome$name
varspl <- selectVar(sgccda.res, block = 'spl', comp = 1)$spl$name
varile <- selectVar(sgccda.res, block = 'ile', comp = 1)$ile$name

#check whether the correlation between components from each data set has been maximized as specified in the design matrix
## For all diets
#col.per.group <- c("darkgreen", "gold", "dodgerblue", "red")
#names(col.per.group) <- levels(Y)
# for ALDR
col.per.group <- c("dodgerblue", "red")
names(col.per.group) <- levels(Y)

pdf(paste0(path,"datasets_correlation.pdf"))
plotDiablo(sgccda.res, ncomp = 1, col.per.group = col.per.group)
dev.off()

pdf(paste0(path,"PCindiv.pdf"))
plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO', ellipse = TRUE, col.per.group = col.per.group)
dev.off()

pdf(paste0(path,"PCcombinedArrows.pdf"))
plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO', col.per.group = col.per.group, ellipse = TRUE)
dev.off()

pdf(paste0(path,"VarContributions.pdf"))
plotLoadings(sgccda.res, comp = 1, contrib = 'max', method = 'median', legend.color = col.per.group, legend = FALSE)
dev.off()

network.res <- network(sgccda.res, save = "pdf", name.save = "myNetwork")
# Save for cytoscape
library(igraph)
write.graph(network.res$gR, file = "myNetwork.gml", format = "gml")

# Calculate performance of the model
set.seed(123)# for reproducibility, only when the `cpus' argument is not used
perf.diablo = perf(sgccda.res, validation = 'Mfold', M = 4, nrepeat = 10, 
                   dist = 'max.dist')
#perf.diablo  # lists the different outputs

# Performance with Majority vote
perf.diablo$MajorityVote.error.rate

pdf(paste0(path,"AUROC.pdf"))
auc.splsda = auroc(sgccda.res, roc.block = "microbiome", roc.comp = 1)
dev.off()

cimDiablo(sgccda.res, save = "pdf", name.save = "cimDiablo", color.blocks= c('lightblue', "darkgreen", 'chocolate3'), color.Y = col.per.group, margins=c(10,5))

pdf(paste0(path,"circos.pdf"))
circosPlot(sgccda.res, cutoff = 0.7, line = TRUE, 
           color.blocks= c('lightblue', "darkgreen", 'chocolate3'),
           color.cor = c("lightgreen","grey20"), size.labels = 1.5, color.Y = col.per.group)
dev.off()

