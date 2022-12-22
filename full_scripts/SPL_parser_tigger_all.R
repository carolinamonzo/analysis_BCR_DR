suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(tigger))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(scoper))
suppressPackageStartupMessages(library(argparser))

# Create a parser
parser <- arg_parser("Script for identifying clone cutoff in IG-Seq spleen")

# Add command line arguments
parser <- add_argument(parser, "--ighv", help = "Path to imgt_mouse_IGHV.fasta")
parser <- add_argument(parser, "--plots_path", help = "Path to place to store plots")
parser <- add_argument(parser, "--to_study", help = "Name of file to store info on samples to further study")
parser <- add_argument(parser, "--log", help = "Name of file for results log")
parser <- add_argument(parser, "--cmd_fof", help = "Path and name of file where to write the CMD for DefineClones")
parser <- add_argument(parser, "--path", help = "Path to files")
parser <- add_argument(parser, "--file", help = "Path and name of file to study")

args <- parse_args(parser)

# Load database for IGHV genes

# EX: "/beegfs/group_lp/home/CMonzo/CM_IGseq/metadata/databases/germlines/imgt/mouse/vdj/imgt_mouse_IGHV.fasta"

ighv <- readIgFasta(args$ighv)

# Set paths
# EX: "/beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/plots/productive_vdj/"
plots_path <- args$plots_path

# "/beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/productive_vdj/"
path <- args$path

# EX: "/beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/productive_vdj/samples_to_study_tigger.log"
repfile <- paste0(path, args$to_study)

# "/beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/productive_vdj/results_tigger.log"
logfile <- paste0(path, args$log)

# "/beegfs/group_lp/home/CMonzo/CM_IGseq/scripts/CMD_FOF/cmd_ILE_defineclones.fof"
pyfile <- args$cmd_fof

# Load files
file <- args$file

## Loop through all samples


# Write in our log file which sample we are running
name <- gsub("_igblast_db-pass_parse-select_parse-select.tsv", "", file)
name <- gsub(path, "", name)

DF <- read_airr(file)
# Remove na in juction length
db <- DF[DF$junction_length != 0,]
# find novel alleles
nv <- tryCatch(
  {
    findNovelAlleles(db, germline_db = ighv, nproc = 50)
  },
  error = function(cond){
    return(findNovelAlleles(db, germline_db = ighv, nproc = 50, germline_min = 20))
  }
)  

# Select novel alleles
if(nrow(selectNovel(nv)) !=0) {
  write("FILE TO LOOK FURTHER INTO!!!", file = repfile, append=TRUE)
  write(name, file=repfile, append=TRUE)
  }

# Now that we have the novel alleles, we genotype and include novel V gene alleles
geno_bayesian <- inferGenotypeBayesian(db, germline_db = ighv, novel = nv, find_unmutated = TRUE)
# Plot the genotype
plotGenotype(geno_bayesian, text_size=10)
ggsave(paste0(plots_path, name, "_Genotype.svg"), plot = last_plot(), device = NULL, path = NULL, width = 20, height = 40, units = "cm")
# save the genotype information in .fasta format to be used later with CreateGermlines.py
gtseq <- genotypeFasta(geno_bayesian, ighv, nv)
writeFasta(gtseq, paste0(path, "results_tigger/", name, "_vgenotype.fasta"))

## Correcting allele calls (reassign gene calls)
# Use the personalized genotype to determine corrected allele assignments
sample_db <- reassignAlleles(db, gtseq)

# Store
write_airr(sample_db, paste0(path, "results_tigger/", name, "_genotype_airr.tsv"))
# Find the set of alleles in the original calls that were not in the genotype
not_in_genotype <- sample_db$v_call %>% strsplit(",") %>% unlist() %>% unique() %>% setdiff(names(gtseq))

## To determine the threshold we use scoper, if it fails, density

res <- tryCatch(
  {
    ## Get clones using spectral clustering from scoper
    threshold <- spectralClones(sample_db, method="vj",
                        germline="germline_alignment")
    # Document in log
    write(name, file=logfile, append=TRUE)
    write("Spectral clustering", file=logfile, append=TRUE)
    # plot it
    plot(threshold, binwidth = 0.02)
    ggsave(paste0(plots_path, name, "_clones_Hist.svg"), plot = last_plot(), device = NULL, path = NULL, width = 20, height = 20, units = "cm")
    # Get threshold value
    slot(threshold, 'eff_threshold')
  },
    error = function(cond){
    # Get clones using hierarchical clustering from shazam
    # Calculate distance
    db_t <- distToNearest(sample_db, model = "ham", normalize = "len", 
                        vCallColumn = "v_call_genotyped", nproc = 50, first = FALSE)
    # Determine Threshold
    threshold <- findThreshold(db_t$dist_nearest, method = "density")
    write(name, file=logfile, append=TRUE)
    write("Hierarchical clustering", file=logfile, append = TRUE)
    # Plot it
    plot(threshold, binwidth = 0.02)
    ggsave(paste0(plots_path, name, "_clones_Hist.svg"), plot = last_plot(), device = NULL, path = NULL, width = 20, height = 20, units = "cm")
    # Get threshold value
    return(round(threshold@threshold, 2))
  }
)

# Now run DefineClones
x <- sprintf('DefineClones.py -d /beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/productive_vdj/results_tigger/%s_genotype_airr.tsv --vf v_call_genotyped --act set --model ham --norm len --dist %s --format airr --nproc 60 --outdir /beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/clonal_groups/', name, res)
  
write(x, file=pyfile, append=TRUE)
