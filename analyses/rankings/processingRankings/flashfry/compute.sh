#!/bin/bash                                            
#SBATCH --mem=20G
i=$SLURM_ARRAY_TASK_ID

module load Java/1.8.0_144
cd /gstore/data/omni/crispr/flashfry


gene=`awk "NR==$i" genes.txt`
database=cas9_human
mkdir outputs/$gene
input=inputs/$gene
input=$input.fasta
output=outputs/$gene
output=$output.output
output_score=$output.scored

java -Xmx4g -jar FlashFry-assembly-1.15.jar \
 discover \
 --database $database \
 --fasta $input \
 --output $output

java -Xmx4g -jar FlashFry-assembly-1.15.jar \
 score \
 --input $output \
 --output $output_score \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot,rank \
 --database $database

