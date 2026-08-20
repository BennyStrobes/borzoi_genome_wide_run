#!/bin/bash
#SBATCH -t 0-0:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB




borzoi_h5_file="${1}"
borzoi_target_file="${2}"
output_tsv_file="${3}"



source ~/.bashrc
conda activate borzoi

python convert_h5_file_to_variant_gene_effects_tsv.py ${borzoi_h5_file} ${borzoi_target_file} ${output_tsv_file}
