#!/bin/bash
#SBATCH -t 0-0:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB 




borzoi_full_h5_file="${1}"
borzoi_gtex_only_h5_file="${2}"
borzoi_gtex_only_target_file="${3}"
chunk_num="${4}"


echo "CHUNK"${chunk_num}

source ~/.bashrc
conda activate borzoi

python subset_h5_file_to_gtex_only.py ${borzoi_full_h5_file} ${borzoi_gtex_only_h5_file} ${borzoi_gtex_only_target_file}