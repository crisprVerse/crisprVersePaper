#!/bin/bash                                            
#SBATCH --mem=10G
i=$SLURM_ARRAY_TASK_ID
cd /gstore/data/omni/crispr/cctop

gene=`awk "NR==$i" genes.txt`
input=outputs/$gene/$gene.xls
output=outputs/$gene/$gene.txt
grep -P "T[1-9]+" $input > $output
