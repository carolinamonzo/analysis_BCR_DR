suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(tigger))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(argparser))


# Load database for IGHV genes
ighv <- readIgFasta("/beegfs/group_lp/home/CMonzo/CM_IGseq/metadata/databases/germlines/imgt/mouse/vdj/imgt_mouse_IGHV.fasta")
# See first germline
# Load files

plots_path <- "/beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/plots/productive_vdj/"

path <- "/beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/productive_vdj/"


## Load sample with Novel alleles. 
print("SAMPLE SGRO-0975, L19251_I8_S8_001_PRIMER-S8_R1 \n")

# Load sample
DF <- read_airr(paste0(path, "dedup-mpg_L19251_I8_S8_001_PRIMER-S8_R1_igblast_db-pass_parse-select_parse-select.tsv"))
# Remove na in juction length
db <- DF[DF$junction_length != 0,]
# find novel alleles
nv <- findNovelAlleles(db, germline_db = ighv, nproc = 20)

# Plot novel alleles:
plotNovel(db, selectNovel(nv)[1,], multiplot = TRUE, v_call = "v_call", j_call = "j_call", seq = "sequence_alignment", junction = "juction", junction_length = "junction_length")
ggsave(paste0(plots_path, "SGRO-0975_IGHV1-5-01_A73G_A87C_A159T_G162C_C249T.svg"), plot = last_plot(), device = NULL, path = NULL, width = 20, height = 20, units = "cm")

plotNovel(db, selectNovel(nv)[2,], multiplot = TRUE, v_call = "v_call", j_call = "j_call", seq = "sequence_alignment", junction = "juction", junction_length = "junction_length")
ggsave(paste0(plots_path, "SGRO-0975_IGHV1-69-01_T226G.svg"), plot = last_plot(), device = NULL, path = NULL, width = 20, height = 20, units = "cm")

#### Genotype and include novel V gene alleles


# Infer the individual's genotype using the bayesian method
geno_bayesian <- inferGenotypeBayesian(db, germline_db = ighv, novel = nv, find_unmutated = TRUE)
# Visualize the begining of the genotype and sequence counts
geno_bayesian %>% arrange(total) %>% slice(1:3)

# In this plots, each row is a gene, with colored cells indicating each of the alleles for that gene that are included in the inferred genotype.

# Plotting the genotypes
plotGenotype(geno_bayesian, text_size=10)
ggsave(paste0(plots_path, "SGRO-0975_Genotype.svg"), plot = last_plot(), device = NULL, path = NULL, width = 20, height = 40, units = "cm")

# save the genotype information in .fasta format to be used later with CreateGermlines.py
gtseq <- genotypeFasta(geno_bayesian, ighv, nv)
writeFasta(gtseq, paste0(path, "results_tigger/mpg_L19251_I8_S8_001_PRIMER-S8_R1_vgenotype.fasta"))


##### Correcting allele calls (reassign gene calls)


# Use the personalized genotype to determine corrected allele assignments
# Updated genotype will be place din the v_call genotyped column
sample_db <- reassignAlleles(db, gtseq)

# Store
write_airr(sample_db, paste0(path, "results_tigger/mpg_L19251_I8_S8_001_PRIMER-S8_R1_genotype_airr.tsv"))
# show some of the corrected gene calls
sample_db %>% filter(v_call != v_call_genotyped) %>%
sample_n(3) %>% select(v_call, v_call_genotyped)


# Find the set of alleles in the original calls that were not in the genotype
not_in_genotype <- sample_db$v_call %>%
  strsplit(",") %>%
  unlist() %>%
  unique() %>%
  setdiff(names(gtseq))

# Determine the fraction of calls that were ambiguous before/after correction and the fraction
# that contained original calls to non-genotype alleles. Note that by dessign, only genotype
# alleles are allowed in the "after" calls

data.frame(Ambiguous = c(mean(grepl(",", sample_db$v_call)),
                         mean(grepl(",", sample_db$v_call_genotyped))),
            NotInGenotype = c(mean(sample_db$v_call %in% not_in_genotype),
                              mean(sample_db$v_call_genotyped %in% not_in_genotype)),
            row.names = c("Before", "After")) %>%
             t() %>% round(3)

evidence <- generateEvidence(sample_db, nv, geno_bayesian, gtseq, ighv, fields = NULL)
evidence %>% select(gene, allele, polymorphism_call, sequences, unmutated_frequency)
