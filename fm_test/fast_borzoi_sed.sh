#!/bin/bash
#SBATCH --gpus 1                             # Request one core
#SBATCH -t 0-10:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-gpu                           # Partition to run in
#SBATCH --mem=22GB  



source ~/.bashrc
conda activate borzoi
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

output_dir="${1}"
vcf_input_file="${2}"
borzoi_training_dir="${3}"
model_num="${4}"


echo $output_dir

python "fast_borzoi_sed.py" -o ${output_dir}"_results_fast.h5" --batch_size 4 --rc --stats logSED,refLog,altLog -t ${borzoi_training_dir}"targets_human.txt" ${borzoi_training_dir}"params_pred.json" ${borzoi_training_dir}"model0_best_f3c"${model_num}".h5" $vcf_input_file
