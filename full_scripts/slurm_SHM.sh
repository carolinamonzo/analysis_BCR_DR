#!/bin/bash

#SBATCH --job-name=SHMS # Job name
#SBATCH --partition=bighead # Partition
#SBATCH --output=./logs/SPL_SHM.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

SING="singularity exec /beegfs/group_lp/home/CMonzo/software/bioinformatics_software_CMonzo.v1.0.sif /bin/bash"

${SING} << EOF
#!/bin/bash
source ~/.bashrc
module load rlang
ulimit -s unlimited
Rscript mutational_load_quantification.R 
EOF

exit
