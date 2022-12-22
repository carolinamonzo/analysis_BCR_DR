#!/bin/bash

SING="singularity exec /beegfs/common/singularity/bioinformatics_software.v3.0.2.sif /bin/bash"


${SING} << EOF
#!/bin/bash
source ~/.bashrc  

module load rlang/3.6.3

Rscript ../full_scripts/SPL_parser_novelAlleles_all.R --ighv /beegfs/group_lp/home/CMonzo/CM_IGseq/metadata/databases/germlines/imgt/mouse/vdj/imgt_mouse_IGHV.fasta --plots_path /beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/plots/productive_vdj/novel_alleles/ --path /beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/productive_vdj/ --file /beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/productive_vdj/dedup-mpg_L19247_S4_S4_001_PRIMER-S4_R1_igblast_db-pass_parse-select_parse-select.tsv --log /beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/productive_vdj/novel_alleles/dedup-mpg_L19247_S4_S4_001_PRIMER-S4_R1_novel_alleles_NEW.log

EOF

exit
