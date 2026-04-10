#####################
# Input data
#####################


# Directory containing pre-trained borzoi models
borzoi_training_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/borzoi_input_data/models/"

# eQTL summary statistics
gtex_v10_eqtl_sumstats_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/eqtl_sumstats/"

# Gtex v10 protein coding genes
gtex_v10_pc_genes_gtf="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v39.gtex.protein_coding.genes.gtf"


# Borzoi target fiel
borzoi_target_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/borzoi_input_data/models/targets_human.txt"

# GTEx sample attributes files
# Contains tissue identity information
gtex_sample_attributes_file="/lab-share/CHIP-Strober-e2/Public/GTEx/gtex_sample_attributes/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"


#####################
# Output data
#####################

# Output root
output_root="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/genome_wide/"

# Directory containing variant-gene pairs and variants to test
variant_gene_dir=${output_root}"input_variant_gene_pairs/"

# Directory containing variant-gene pairs and variants to test
borzoi_pred_dir=${output_root}"borzoi_predictions/"



variant_gene_pair_file=${variant_gene_dir}"variant_gene_pairs_to_test.txt"
variant_output_stem=${variant_gene_dir}"variants_to_test_"
if false; then
sh extract_variant_gene_pairs_to_test.sh $gtex_v10_eqtl_sumstats_dir $variant_gene_pair_file $variant_output_stem $gtex_v10_pc_genes_gtf
fi


borzoi_gtex_only_target_file=${borzoi_pred_dir}"targets_gtex_only_ordered.txt"
if false; then
source ~/.bashrc
conda activate borzoi
python process_borzoi_target_files_for_gtex_only_targets.py $borzoi_target_file $borzoi_gtex_only_target_file $gtex_sample_attributes_file
fi

model_num="0"
if false; then
for chunk_num in {0..29}
do
	variant_vcf_file=$variant_output_stem"chunked_variants_"${chunk_num}".vcf"
	sbatch fast_borzoi_sed.sh $borzoi_pred_dir"model_"${model_num}"_chunk_"${chunk_num}"_borzoi_results.h5" ${variant_vcf_file} $borzoi_training_dir $model_num $variant_gene_pair_file
done
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




