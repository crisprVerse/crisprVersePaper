cd ..
mkdir tmp
ml Java/1.8.0_144 
java -Xmx4g -jar FlashFry-assembly-1.15.jar \
 index \
 --tmpLocation ./tmp \
 --database cas9_human \
 --reference /gstore/data/omni/crispr/crisprIndices/genomes/hg38/hg38.fa \
 --enzyme spcas9ngg


mkdir tmp
 java -Xmx4g -jar FlashFry-assembly-1.15.jar \
 index \
 --tmpLocation ./tmp \
 --database cas9_human \
 --reference /Users/fortinj2/crisprIndices/genomes/hg38/hg38.fa \
 --enzyme spcas9ngg
