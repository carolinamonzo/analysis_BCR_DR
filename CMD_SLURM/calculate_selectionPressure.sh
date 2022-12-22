#!/bin/bash

#SBATCH --job-name=SAgeSelection # Job name
#SBATCH --partition=hooli # Partition # Run on one node
#SBATCH --output=./logs/SPL_AGE_selectionPressure.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

SING="singularity exec /beegfs/common/singularity/bioinformatics_software.v3.0.2.sif /bin/bash"

${SING} << EOF
#!/bin/bash
source ~/.bashrc
module load rlang
ulimit -s unlimited
Rscript ../full_scripts/Selection_pressure.R 
EOF

exit
