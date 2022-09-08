#!/bin/bash           
#SBATCH -N 1                                 
#SBATCH --qos=short                       
#SBATCH --mem=10G
#SBATCH -J extract
#SBATCH -o extract_%A_%a.out                
#SBATCH -e extract_%A_%a.err                
#SBATCH --chdir=/gstore/data/omni/crispr/cctop
#SBATCH --array 1-20395

./extract.sh