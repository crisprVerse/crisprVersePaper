chmod -R 700 compute*
sbatch compute_array.sh 

chmod -R 700 extract*
sbatch extract_array.sh 
