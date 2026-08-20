#!/bin/bash
#SBATCH -t 0-0:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB 




borzoi_full_h5_file="${1}"
borzoi_specified_indices_h5_file="${2}"
borzoi_specified_indices_target_file="${3}"



source ~/.bashrc
conda activate borzoi

python subset_h5_file_to_specified_indices.py ${borzoi_full_h5_file} ${borzoi_specified_indices_h5_file} ${borzoi_specified_indices_target_file}