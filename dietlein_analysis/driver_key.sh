#####################
# Input data
#####################


# Directory containing pre-trained borzoi models
borzoi_training_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/borzoi_input_data/models/"


# Gtex v10 protein coding genes
gtex_v10_pc_genes_gtf="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v39.gtex.protein_coding.genes.gtf"

# Borzoi target file
borzoi_gtex_only_target_file="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/genome_wide/borzoi_predictions/targets_gtex_only_ordered.txt"

# Variant vcf file
variant_vcf_file="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/dietlein_analysis/input_variant_gene_pairs/sarc_col31_test_hg38.vcf"


#####################
# Output data
#####################

# Output root
output_root="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/dietlein_analysis/"

# Directory containing variant-gene pairs and variants to test
variant_gene_dir=${output_root}"input_variant_gene_pairs/"

# Directory containing variant-gene pairs and variants to test
borzoi_pred_dir=${output_root}"borzoi_predictions/"




model_num="0"
if false; then
sbatch fast_borzoi_sed.sh $borzoi_pred_dir"model_"${model_num}"dietlein_analysis_borzoi_results.h5" ${variant_vcf_file} $borzoi_training_dir $model_num
fi







model_num="0"
if false; then
for chunk_num in {0..29}
do
	borzoi_full_h5_file=$borzoi_pred_dir"model_"${model_num}"_chunk_"${chunk_num}"_borzoi_results.h5"
	borzoi_gtex_only_h5_file=$borzoi_pred_dir"model_"${model_num}"_chunk_"${chunk_num}"_borzoi_gtex_only_results.h5"
	sbatch subset_h5_file_to_gtex_only.sh ${borzoi_full_h5_file} ${borzoi_gtex_only_h5_file} ${borzoi_gtex_only_target_file} ${chunk_num}
done
fi



model_num="0"
if false; then
for chunk_num in {0..29}
do
	borzoi_full_h5_file=$borzoi_pred_dir"model_"${model_num}"_chunk_"${chunk_num}"_borzoi_results.h5"
	borzoi_non_gtex_h5_file=$borzoi_pred_dir"model_"${model_num}"_chunk_"${chunk_num}"_borzoi_non_gtex_interesting_results.h5"
	sbatch subset_h5_file_to_gtex_only.sh ${borzoi_full_h5_file} ${borzoi_non_gtex_h5_file} ${borzoi_non_gtex_target_file} ${chunk_num}
done
fi













