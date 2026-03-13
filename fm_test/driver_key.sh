#!/bin/bash
#SBATCH -t 0-1:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)

#################
# Input data
#################
gtex_target_file="/lab-share/CHIP-Strober-e2/Public/Borzoi_data/w5/gtex_targets.txt"

baskerville_code_dir="/home/ch271704/tools/baskerville/src/baskerville/scripts/"

borzoi_micro_json_file="/home/ch271704/tools/borzoi/tutorials/latest/train_model/params_gtex_micro.json"
borzoi_micro_json_file="/home/ch271704/tools/borzoi/tutorials/latest/train_model/params_gtex_micro_100_iter.json"

raw_fine_mapping_eqtl_results_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/fine_mapping/v10/"

gtex_sample_attributes_file="/lab-share/CHIP-Strober-e2/Public/GTEx/gtex_sample_attributes/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"

eqtl_sumstats_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/eqtl_sumstats/"

borzoi_code_dir="/home/ch271704/tools/borzoi/src/scripts/"

protein_coding_genes_file="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gene_tss.bed"

gene_tss_file="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/genecode.v39.GRCh38.bed"

gtex_v10_gene_gtf_file="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v39.GRCh38.genes.gtf"

borzoi_training_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/borzoi_input_data/models/"

raw_fine_mapping_file="/lab-share/CHIP-Strober-e2/Public/GTEx/fine_mapping/v8/GTEx_49tissues_release1.tsv"

gtex_v8_expression_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/GTEx_Analysis_v8_eQTL_expression_matrices/"


#################
# Do fine-mapping analysis first!
#################
# Output data
#################
output_root="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/fine_mapped_eqtls/"

processed_fm_results_dir=${output_root}"processed_fine_mapped_data/"

borzoi_fm_pred_default_dir=${output_root}"borzoi_fm_pred_default_code/"



#################
# Run code
#################
# Re-organize fine-mapped eQTL results
cross_tissue_fine_mapping_summary_file=${processed_fm_results_dir}"fine_mapping_summary.txt"
if false; then
source ~/.bashrc
conda activate borzoi
python organize_fine_mapped_eqtls.py $raw_fine_mapping_file $cross_tissue_fine_mapping_summary_file $gtex_v8_expression_dir
fi


# convert to vcf
fm_variant_vcf_file=${processed_fm_results_dir}"fine_mapping_variant_vcf.txt"
if false; then
source ~/.bashrc
conda activate borzoi
python convert_fm_eqtls_to_variant_vcf.py $cross_tissue_fine_mapping_summary_file $fm_variant_vcf_file
fi
small_fm_variant_vcf_file=${processed_fm_results_dir}"fine_mapping_variant_vcf_small.txt"


model_num="0"
if false; then
sbatch borzoi_sed.sh $borzoi_fm_pred_default_dir"model_"${model_num}"_"${split_number} ${small_fm_variant_vcf_file} $borzoi_training_dir $model_num
fi


model_num="0"
if false; then
sh fast_borzoi_sed.sh $borzoi_fm_pred_default_dir"model_"${model_num} ${small_fm_variant_vcf_file} $borzoi_training_dir $model_num
fi



echo "CONSIDER RUNNING WITH GTEX fasta!"


if false; then
source ~/.bashrc
conda activate borzoi
python compare_results.py $borzoi_fm_pred_default_dir"model_"${model_num}
fi
