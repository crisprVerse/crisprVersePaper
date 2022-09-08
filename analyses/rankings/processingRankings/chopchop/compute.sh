#!/bin/bash                                            
#SBATCH --mem=20G
i=$SLURM_ARRAY_TASK_ID

cd /gstore/data/omni/crispr/chopchop
source /gstore/data/omni/crispr/python/miniconda3/etc/profile.d/conda.sh
conda activate chopchop

gene=`awk "NR==$i" genes.txt`
input=../flashfry/inputs/$gene.fasta
outdir=outputs/$gene
outfile=$outdir/results.txt
mkdir $outdir
./chopchop/chopchop.py -F -Target $input -G hg38 --scoringMethod DOENCH_2016  -o $outdir > $outfile
cd ./$outdir
rm -rf *.offtargets
rm -rf sequence.fa
cd ../..
