#####################
# Input data
#####################


# Directory containing pre-trained borzoi models
borzoi_training_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/borzoi_input_data/models/"

# eQTL summary statistics
gtex_v10_eqtl_sumstats_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/eqtl_sumstats/"

# Gtex v10 protein coding genes
gtex_v10_pc_genes_gtf="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v39.gtex.protein_coding.genes.gtf"


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

model_num="0"
chunk_num="0"
variant_vcf_file=$variant_output_stem"chunked_variants_"${chunk_num}".vcf"
if false; then
sbatch fast_borzoi_sed.sh $borzoi_pred_dir"model_"${model_num}"_chunk_"${chunk_num}"_borzoi_results.h5" ${variant_vcf_file} $borzoi_training_dir $model_num $variant_gene_pair_file
fi
