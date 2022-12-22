#!/bin/bash

#SBATCH --job-name=ILEisotypes # Job name
#SBATCH --partition=himem # Partition
#SBATCH --nodes=1 # Run on one node
#SBATCH --output=./ILEalpha_isotypes.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

SING="singularity exec /beegfs/common/singularity/bioinformatics_software.v3.0.2.sif /bin/bash"

${SING} << EOF
#!/bin/bash
source ~/.bashrc

module load rlang/3.6.3

ulimit -s unlimited
Rscript isotypes_alpha_calculate.R

EOF

exit
