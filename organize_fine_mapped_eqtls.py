import numpy as np
import os
import sys
import pdb
import gzip

def create_dictionary_list_of_analyzed_gene_tissue_pairs(gtex_expression_dir):
	dicti = {}
	for file_name in os.listdir(gtex_expression_dir):
		if file_name.endswith('expression.bed.gz') == False:
			continue

		full_file_name = gtex_expression_dir + file_name
		tissue_name = file_name.split('.v8')[0]

		f = gzip.open(full_file_name, 'rt')
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			gene_id1 = data[3]
			gene_id2 = data[3].split('.')[0]

			dicti[gene_id1 + ':' + tissue_name] = 1
			dicti[gene_id2 + ':' + tissue_name] = 1
		f.close()
	return dicti

def get_ordered_tissue_names(raw_fine_mapping_file):
	f = open(raw_fine_mapping_file)
	head_count = 0
	tissues = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissues[data[10]] = 1
	f.close()

	ordered_tissues = np.sort([*tissues])


	return ordered_tissues


def extract_variant_gene_pair_info(raw_fine_mapping_file, tissue_name_to_index, gene_tissue_list,ordered_tissue_names):
	dicti = {}
	n_tiss = len(tissue_name_to_index)

	f = open(raw_fine_mapping_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count==0:
			head_count = head_count + 1
			continue
		if data[9] != 'SUSIE':
			continue
		if data[8] != 'GTEx':
			continue
		tissue_name = data[10]
		variant_id = data[4]
		gene_id = data[11]
		posterior_effect = float(data[-2])

		vg_pair = variant_id + ':' + gene_id
		if vg_pair not in dicti:
			dicti[vg_pair] = np.asarray([np.nan]*n_tiss)
		dicti[vg_pair][tissue_name_to_index[tissue_name]] = posterior_effect
	f.close()

	vg_pairs = [*dicti]

	for vg_pair in vg_pairs:
		gene_id = vg_pair.split(':')[1]
		for tiss_iter, tissue_name in enumerate(ordered_tissue_names):
			if np.isnan(dicti[vg_pair][tiss_iter]):
				gt = gene_id + ':' + tissue_name
				if gt in gene_tissue_list:
					dicti[vg_pair][tiss_iter] = 0.0
			else:
				gt = gene_id + ':' + tissue_name
				if gt not in gene_tissue_list:
					print('assumption error')
					pdb.set_trace()

	return dicti

######################
# Command line args
######################
raw_fine_mapping_file = sys.argv[1]
cross_tissue_fine_mapping_summary_file = sys.argv[2] # Output file
gtex_expression_dir = sys.argv[3]


####################
# First create dictionary list of gene-tissue pairs tested
gene_tissue_list = create_dictionary_list_of_analyzed_gene_tissue_pairs(gtex_expression_dir)

####################
# Second get ordered list of tissue names
ordered_tissue_names = get_ordered_tissue_names(raw_fine_mapping_file)
# Create mapping from tissue name to index
tissue_name_to_index = {}
for tiss_iter, tissue_name in enumerate(ordered_tissue_names):
	tissue_name_to_index[tissue_name] = tiss_iter


####################
# Third extract variant_gene_pair info
vg_pair_info = extract_variant_gene_pair_info(raw_fine_mapping_file, tissue_name_to_index, gene_tissue_list,ordered_tissue_names)


####################
# print to output
t = open(cross_tissue_fine_mapping_summary_file,'w')
# Header
t.write('variant_id\tgene_id\t' + '\t'.join(ordered_tissue_names) + '\n')

for vg_pair in np.sort([*vg_pair_info]):
	variant_id, gene_id = vg_pair.split(':')

	vec = vg_pair_info[vg_pair]

	t.write(variant_id + '\t' + gene_id + '\t' + '\t'.join(vec.astype(str)) + '\n')

t.close()
