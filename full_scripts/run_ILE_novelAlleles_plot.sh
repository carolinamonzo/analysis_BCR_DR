#!/bin/bash

SING="singularity exec /beegfs/common/singularity/bioinformatics_software.v3.0.2.sif /bin/bash"


${SING} << EOF
#!/bin/bash
source ~/.bashrc  

module load rlang/3.6.3

Rscript ILE_novelAlleles_plot.R

EOF

exit
