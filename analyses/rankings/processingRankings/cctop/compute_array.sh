#!/bin/bash           
#SBATCH -N 1                                 
#SBATCH --qos=medium                       
#SBATCH --mem=20G
#SBATCH -J compute
#SBATCH -o compute_%A_%a.out                
#SBATCH -e compute_%A_%a.err                
#SBATCH --chdir=/gstore/data/omni/crispr/cctop
#SBATCH --array 1-20395

./compute.sh