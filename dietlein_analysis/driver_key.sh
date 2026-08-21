#####################
# Input data
#####################


# Directory containing pre-trained borzoi models
borzoi_training_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/borzoi_input_data/models/"


# Gtex v10 protein coding genes
gtex_v10_pc_genes_gtf="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v39.gtex.protein_coding.genes.gtf"

# Borzoi target file
borzoi_gtex_only_target_file="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/genome_wide/borzoi_predictions/targets_gtex_only_ordered.txt"
borzoi_fibroblast_only_target_file="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/genome_wide/borzoi_predictions/targets_fibroblast_only_ordered.txt"
borzoi_dietlein_only_target_file="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/genome_wide/borzoi_predictions/targets_dietlein_cell_types_only_ordered.txt"
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
borzoi_full_h5_file=$borzoi_pred_dir"model_"${model_num}"dietlein_analysis_borzoi_results.h5"
borzoi_fibroblast_only_h5_file=$borzoi_pred_dir"model_"${model_num}"dietlein_analysis_borzoi_fibroblast_only_results.h5"
if false; then
sh subset_h5_file_to_specified_indices.sh ${borzoi_full_h5_file} ${borzoi_fibroblast_only_h5_file} ${borzoi_fibroblast_only_target_file}
fi

if false; then
# Make tsv file where each row is a variant-gene pair along with predicted effects
borzoi_fibroblast_only_effects_tsv_file=$borzoi_pred_dir"model_"${model_num}"dietlein_analysis_borzoi_fibroblast_only_variant_gene_effects.tsv"
sh convert_h5_file_to_variant_gene_effects_tsv.sh ${borzoi_fibroblast_only_h5_file} ${borzoi_fibroblast_only_target_file} ${borzoi_fibroblast_only_effects_tsv_file}

# Filter tsv to variant-gene pairs for COL1A2 (ENSG00000164692) and COL3A1 (ENSG00000168542)
borzoi_fibroblast_only_COL1A2_COL3A1_effects_tsv_file=$borzoi_pred_dir"model_"${model_num}"dietlein_analysis_borzoi_fibroblast_only_COL1A2_COL3A1_only_variant_gene_effects.tsv"
python filter_variant_gene_effects_tsv_to_specified_genes.py ${borzoi_fibroblast_only_effects_tsv_file} ${borzoi_fibroblast_only_COL1A2_COL3A1_effects_tsv_file} "ENSG00000164692.19,ENSG00000168542.16"
fi


if false; then
model_num="0"
borzoi_full_h5_file=$borzoi_pred_dir"model_"${model_num}"dietlein_analysis_borzoi_results.h5"
borzoi_other_cell_types_only_h5_file=$borzoi_pred_dir"model_"${model_num}"dietlein_analysis_borzoi_other_cell_types_only_results.h5"
sh subset_h5_file_to_specified_indices.sh ${borzoi_full_h5_file} ${borzoi_other_cell_types_only_h5_file} ${borzoi_dietlein_only_target_file}



# Make tsv file where each row is a variant-gene pair along with predicted effects
borzoi_other_cell_types_only_effects_tsv_file=$borzoi_pred_dir"model_"${model_num}"dietlein_analysis_borzoi_other_cell_types_only_variant_gene_effects.tsv"
sh convert_h5_file_to_variant_gene_effects_tsv.sh ${borzoi_other_cell_types_only_h5_file} ${borzoi_dietlein_only_target_file} ${borzoi_other_cell_types_only_effects_tsv_file}

# Filter tsv to variant-gene pairs for COL1A2 (ENSG00000164692) and COL3A1 (ENSG00000168542)
borzoi_other_cell_types_only_COL1A2_COL3A1_effects_tsv_file=$borzoi_pred_dir"model_"${model_num}"dietlein_analysis_borzoi_other_cell_types_only_COL1A2_COL3A1_only_variant_gene_effects.tsv"
python filter_variant_gene_effects_tsv_to_specified_genes.py ${borzoi_other_cell_types_only_effects_tsv_file} ${borzoi_other_cell_types_only_COL1A2_COL3A1_effects_tsv_file} "ENSG00000164692.19,ENSG00000168542.16"
fi
















