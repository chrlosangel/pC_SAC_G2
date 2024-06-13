#!/bin/bash

#SBATCH --job-name=generate_data
#SBATCH --time=10:00:00
#SBATCH --partition=pe2
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=cangel@nygenome.org

# define an array of $2 values
#change if needed
#resolutions=(5000)
#lrf=(4) 
main_dir="/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning"
resolutions=$6
lrf=$7
mkdir -p $main_dir/PBS
# loop over the array of $2 values
for res in "${resolutions[@]}"; do
    for LR in "${lrf[@]}";do

        # specify the output log file
        log_file="$main_dir/PBS/generate_data_${res}_${LR}.log"

        # run the script with the current $2 value 
        # change script if needed 
        srun --ntasks 1 -u python $main_dir/generate_data.py --hic-file $1 \
                            --hic-res $res --chromosome $2 \
                            --start $3  --end $4 \
                            --d-factor $5 --bool-res $8 \
                            --l_res $LR --split_chains $9 > "$log_file" 2>&1
                            
        echo $'\n' >> "$log_file"
        echo "Reading file $1" >> "$log_file"
        echo "Resolution map $res" >> "$log_file"
        echo "Chromosome: $2" >> "$log_file"
        echo "Start coordinate: $3" >> "$log_file"
        echo "End coordinate: $4" >> "$log_file"
        echo "Downsampling factor: $5" >> "$log_file"
        echo "Bool resolution : $8" >> "$log_file"
        echo "Split Chains : $9" >> "$log_file"
        echo "LR map factor (default 4) : $LR" >> "$log_file"
        done
done
