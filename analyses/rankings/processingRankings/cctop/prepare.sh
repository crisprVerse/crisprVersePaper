cp ../flashfry/genes.txt ./
mkdir outputs



cctop --input example.fasta --index /Users/fortinj2/crisprIndices/bowtie/hg38/hg38 --output ./temp
ml bowtie
conda create -n cctop
conda activate cctop
conda install bowtie --channel Bioconda
pip install cctop
export PATH=/gstore/data/omni/crispr/programs/bowtie-1.3.1-src:$PATH
