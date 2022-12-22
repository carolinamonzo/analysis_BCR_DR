#!/bin/bash

#SBATCH --job-name=CLONES
#SBATCH --partition=himem
#SBATCH --output=./logs/ILE_tigger_clones.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cmonzo@age.mpg.de
#SBATCH --cpus-per-task=60

SING="singularity exec /beegfs/common/singularity/bioinformatics_software.v3.0.2.sif /bin/bash"


${SING} << EOF
#!/bin/bash
source ~/.bashrc  

module load rlang/3.6.3

Rscript ILE_tigger_all.R

EOF

exit
