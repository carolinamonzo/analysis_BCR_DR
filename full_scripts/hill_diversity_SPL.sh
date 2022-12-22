#!/bin/bash

#SBATCH --job-name=SPL-CLONES # Job name
#SBATCH --partition=hugemem # Partition
#SBATCH --nodes=1 # Run on one node
#SBATCH --output=./logs/SPL_alpha.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

# Script from Jorge - 23062020

SING="singularity exec /beegfs/common/singularity/bioinformatics_software.v3.0.2.sif /bin/bash"


${SING} << EOF
#!/bin/bash
source ~/.bashrc  
module load rlang/3.6.3
ulimit -s unlimited

R --slave -e 'Cstack_info()["size"]'

Rscript -e '
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(tigger))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(scoper))

names(IG_COLORS) <- c( "IGM", "IGG12", "IGA", "IGD", "IGE", "IGG3")

db <- read_airr("../../analysis/created_germlines/reannotated_db/merged_changeo_SPL.tsv")

sample_colors <- c("5AL_lifelong"="#D3E0F1", "16AL_lifelong"="#A2C2DD", "20AL_lifelong"="#6897C6", "24AL_lifelong"="#3D64A8", "5DR_lifelong"="#ECBAA1", "20DR_lifelong"="#DD694D", "24DR_lifelong"="#AD1E22","20AL_DR16M"="#FFFFA2", "24AL_DR16M"="#EDD853", "24AL_DR20M"="#779D45")

curve <- estimateAbundance(db, group = "biogroup", ci = 0.95, nboot = 200, clone = "clone_id")

plot(curve, colors = sample_colors, legend_title="Sample")
ggsave("/beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/plots/Clonal_abundance_SPL.svg", plot = last_plot(), device = NULL, path=NULL, width = 20, height = 10, units = "cm")

sample_curve <- alphaDiversity(db, group = "biogroup", clone="clone_id", min_q=0, max_q=4, step_q=0.1, ci=0.95, nboot=200)

sample_main <- paste0("Sample diversity")

plot(sample_curve, colors=sample_colors, main_title=sample_main, legend_title="Sample")
ggsave("../../analysis/plots/Alpha_Hill_SPL.svg", plot = last_plot(), device = NULL, path=NULL, width = 20, height = 10, units = "cm")

sample_test <- alphaDiversity(db, group="biogroup", clone="clone_id", min_q=0, max_q=2, step_q=1, nboot=200)

write.table(sample_test@tests, file = "../../analysis/plots/alpha_values_stats_SPL.tsv", append=TRUE, col.names = TRUE, row.names = TRUE)


'

EOF

exit

