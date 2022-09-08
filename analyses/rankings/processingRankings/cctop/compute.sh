#!/bin/bash                                            
#SBATCH --mem=20G
i=$SLURM_ARRAY_TASK_ID
cd /gstore/data/omni/crispr/cctop
conda activate cctop

gene=`awk "NR==$i" genes.txt`
index=/gstore/data/omni/crispr/crisprIndices/bowtie/hg38/hg38
mkdir outputs/$gene
input=../flashfry/inputs/$gene
input=$input.fasta
output=outputs/$gene

cctop --input $input --index $index --output $output