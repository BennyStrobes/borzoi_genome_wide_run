
#####################
# Input data
#####################

# Directory containing pre-trained borzoi models
borzoi_training_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/borzoi_input_data/models/"


# Gtex v10 protein coding genes
gtex_v10_pc_genes_gtf="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v39.gtex.protein_coding.genes.gtf"


# Borzoi target fiel
borzoi_target_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/borzoi_input_data/models/targets_human.txt"
borzoi_gtex_target_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/borzoi_input_data/models/targets_gtex.txt"


# Example variant-gene pairs
# Output root
example_variant_gene_dir="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/genome_wide/input_variant_gene_pairs/"
example_variant_gene_pair_file=${example_variant_gene_dir}"variant_gene_pairs_to_test.txt"
existing_variant_vcf_file=$example_variant_gene_dir"variants_to_test_chunked_variants_0.vcf"


#####################
# Output data
#####################

# Output root
output_root="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/genome_wide/justin_example/"




###################
# Analysis
###################

#####################
# Make small, example vcf containing 500 variants to run Borzoi for
example_variant_vcf_file=$output_root"example_variants.vcf"
head -1000 $existing_variant_vcf_file > $example_variant_vcf_file

######################
# Display first 10 lines of vcf
head -10 $example_variant_vcf_file



######################
# Run Borzoi on these 10 variants

# Load in borzoi package using conda (Justin: you will need to first install your own borzoi package)
source ~/.bashrc
conda activate borzoi
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"


# Location/name of file to write results to 
output_file=${output_root}"example_output2.h5"

# Run Borzoi on variants in $example_variant_vcf_file, while further limiting to variant-gene pairs found in ${example_variant_gene_pair_file}
python "fast_borzoi_sed.py" -o ${output_file} -v ${example_variant_gene_pair_file} --batch_size 4 --rc --stats logSED,refLog,altLog -t ${borzoi_gtex_target_file} ${borzoi_training_dir}"params_pred.json" ${borzoi_training_dir}"model0_best_f3c0.h5" $example_variant_vcf_file












