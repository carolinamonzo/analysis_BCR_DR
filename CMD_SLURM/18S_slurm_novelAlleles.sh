#!/bin/bash

#SBATCH --job-name=CLONES # Job name
#SBATCH --partition=himem # Partition
#SBATCH --nodes=1 # Run on one node
#SBATCH --output=./logs/SPL_novelAlleles.log # Standard output and error log
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=cmonzo@age.mpg.de

SING="singularity exec /beegfs/common/singularity/bioinformatics_software.v3.0.2.sif /bin/bash"

${SING} << EOF
#!/bin/bash
source ~/.bashrc
module load rlang/3.6.3

parallel -j 1 :::: ../CMD_FOF/cmd_SPL_novelAlleles.fof

EOF

exit
