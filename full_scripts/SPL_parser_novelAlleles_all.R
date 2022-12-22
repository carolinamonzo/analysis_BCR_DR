suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(tigger))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(argparser))

# Create a parser
parser <- arg_parser("Script for identifying novel alleles")

# Add command line arguments
parser <- add_argument(parser, "--ighv", help = "Path to imgt_mouse_IGHV.fasta") 
parser <- add_argument(parser, "--plots_path", help = "Path to place to store plots") 
parser <- add_argument(parser, "--path", help = "Path to files")
parser <- add_argument(parser, "--file", help = "Path and name of file to study")
parser <- add_argument(parser, "--log", help = "Name of file for results log")
args <- parse_args(parser)


# Load database for IGHV genes
# EX: ighv <- readIgFasta("/beegfs/group_lp/home/CMonzo/CM_IGseq/metadata/databases/germlines/imgt/mouse/vdj/imgt_mouse_IGHV.fasta")

ighv <- readIgFasta(args$ighv)

# Load files
# EX: "/beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/plots/productive_vdj/"
plots_path <- args$plots_path

# EX: "/beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/productive_vdj/"
path <- args$path

# EX: "/beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/productive_vdj/results_novel_alleles.log"
logfile <- args$log

file <- args$file

# Get name
name <- gsub("_igblast_db-pass_parse-select_parse-select.tsv", "", file)
name <- gsub(path, "", name)

# Load sample
DF <- read_airr(file)
# Remove na in juction length
db <- DF[DF$junction_length != 0,]
# find novel alleles
nv <- findNovelAlleles(db, germline_db = ighv, nproc = 20)

# Get novel alleles
nov <- selectNovel(nv)

write("\nNovel alleles\n", file = logfile, append = TRUE)
write.table(nov, file = logfile, append = TRUE, col.names = TRUE, row.names = TRUE)

# Plot novel alleles:

for (i in 1:dim(nov)[1]) {
    plotNovel(db, nov[i,], multiplot = TRUE, v_call = "v_call", j_call = "j_call", seq = "sequence_alignment", junction = "junction", junction_length = "junction_length")
    ggsave(paste0(plots_path, name, "_allele_", i, ".svg"), plot = last_plot(), device = NULL, path = NULL, width = 20, height = 20, units = "cm")
}

#### Genotype and include novel V gene alleles


# Infer the individual's genotype using the bayesian method
geno_bayesian <- inferGenotypeBayesian(db, germline_db = ighv, novel = nv, find_unmutated = TRUE)

write("\nVisualize the begining of the genotype and sequence counts\n", file = logfile, append = TRUE)
write.table(geno_bayesian %>% arrange(total), file = logfile, append = TRUE, col.names = TRUE, row.names = TRUE)


# In this plots, each row is a gene, with colored cells indicating each of the alleles for that gene that are included in the inferred genotype.

# Plotting the genotypes
plotGenotype(geno_bayesian, text_size=10)
ggsave(paste0(plots_path, name, "_Genotype.svg"), plot = last_plot(), device = NULL, path = NULL, width = 20, height = 40, units = "cm")

# save the genotype information in .fasta format to be used later with CreateGermlines.py
gtseq <- genotypeFasta(geno_bayesian, ighv, nv)
writeFasta(gtseq, paste0(path, "results_tigger/", name, "_vgenotype.fasta"))


##### Correcting allele calls (reassign gene calls)


# Use the personalized genotype to determine corrected allele assignments
# Updated genotype will be place din the v_call genotyped column
sample_db <- reassignAlleles(db, gtseq)

# Store
write_airr(sample_db, paste0(path, "results_tigger/", name, "_genotype_airr.tsv"))
# show some of the corrected gene calls

write("\nShow some of the corrected gene calls\n", file = logfile, append = TRUE)
write.table(sample_db %>% filter(v_call != v_call_genotyped) %>% sample_n(3) %>% select(v_call, v_call_genotyped), file = logfile, append = TRUE, col.names = TRUE, row.names = TRUE)

# Find the set of alleles in the original calls that were not in the genotype
not_in_genotype <- sample_db$v_call %>%
  strsplit(",") %>%
  unlist() %>%
  unique() %>%
  setdiff(names(gtseq))

# Determine the fraction of calls that were ambiguous before/after correction and the fraction
# that contained original calls to non-genotype alleles. Note that by dessign, only genotype
# alleles are allowed in the "after" calls
write("\nShow percentage of ambiguous calls and reduction after reassignment\n", file = logfile, append = TRUE)

write.table(data.frame(Ambiguous = c(mean(grepl(",", sample_db$v_call)),
                         mean(grepl(",", sample_db$v_call_genotyped))),
            NotInGenotype = c(mean(sample_db$v_call %in% not_in_genotype),
                              mean(sample_db$v_call_genotyped %in% not_in_genotype)),
            row.names = c("Before", "After")) %>%
             t() %>% round(3), file = logfile, append = TRUE, row.names = TRUE, col.names = TRUE)

evidence <- generateEvidence(sample_db, nv, geno_bayesian, gtseq, ighv, fields = NULL)
write("\nEvidence, number of sequences unambiguously reassigned\n", file = logfile, append = TRUE)

write.table(evidence %>% select(gene, allele, polymorphism_call, sequences, unmutated_frequency), file = logfile, append = TRUE, col.names = TRUE, row.names = TRUE)
