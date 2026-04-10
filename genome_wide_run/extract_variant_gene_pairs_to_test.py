import numpy as np
import os
import sys
import pdb
import pyarrow.parquet as pq
import pandas as pd
import pyarrow.compute as pc



def extract_ordered_list_of_all_variant_ids_on_this_chromosome(gtex_sumstats_dir, chrom_string, pc_genes):
	vg_pairs = {}
	pc_genes_set = set(pc_genes)

	for file_name in os.listdir(gtex_sumstats_dir):
		if not file_name.endswith('chr' + str(chrom_string) + '.parquet'):
			continue

		print(file_name)
		pf = pq.ParquetFile(gtex_sumstats_dir + file_name)

		for rg in range(pf.num_row_groups):
			table = pf.read_row_group(
				rg,
				columns=['gene_id', 'variant_id', 'tss_distance', 'af']
			)

			if table.num_columns == 0:
				continue

			gene_col = table['gene_id']
			variant_col = table['variant_id']
			dist_col = table['tss_distance']
			af_col = table['af']

			mask = pc.less_equal(pc.abs(dist_col), 100000)

			filt = table.filter(mask)

			if filt.num_rows == 0:
				continue

			gene_ids = filt['gene_id'].to_pylist()
			variant_ids = filt['variant_id'].to_pylist()

			for gene_id, variant_id in zip(gene_ids, variant_ids):
				if gene_id.split('.', 1)[0] not in pc_genes_set:
					continue
				vg_pairs[f"{variant_id}:{gene_id}"] = 1

	return vg_pairs

def extract_dictionary_list_of_protein_coding_genes(pc_genes_gtf):
	f = open(pc_genes_gtf)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		ens_id = data[8].split(';')[0].split('"')[1]
		if ens_id.startswith('ENSG') == False:
			print('assumption oernroro')
			pdb.set_trace()
		dicti[ens_id.split('.')[0]] = 1

	f.close()

	return dicti



def make_variant_vcf_file(variant_gene_pair_file, vcf_file):
	f = open(variant_gene_pair_file)
	head_count = 0
	var_dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		var_dicti[data[0]] = 1
	f.close()

	tupler = []
	for var_id in [*var_dicti]:
		if len(var_id.split('_')) == 1:
			print(var_id)
			continue
		chromer = int(var_id.split('_')[0].split('hr')[1])
		position = int(var_id.split('_')[1])
		tupler.append((chromer, position, var_id))

	tupler.sort(key=lambda x: (x[0], x[1]))

	t = open(vcf_file,'w')

	for tup in tupler:
		var_id = tup[2]
		chromer, pos, a1, a2, garbage = var_id.split('_')

		if len(a1) != 1 or len(a2) != 1:
			continue

		t.write(chromer + '\t' + pos + '\t' + var_id + '\t' + a1 + '\t' + a2 + '\n')

	t.close()

	return





# Command line args
eqtl_sumstats_dir = sys.argv[1]
variant_gene_pair_file = sys.argv[2]
variant_output_stem = sys.argv[3]
pc_genes_gtf = sys.argv[4]

# Extract dictionary list of protein coding genes
pc_genes = extract_dictionary_list_of_protein_coding_genes(pc_genes_gtf)

t = open(variant_gene_pair_file,'w')
t.write('variant_id\tgene_id\n')

for chrom_num in range(1,23):
	print(chrom_num)

	# Extract ordered list of variant ids in this chromosome
	variant_gene_pairs = extract_ordered_list_of_all_variant_ids_on_this_chromosome(eqtl_sumstats_dir, str(chrom_num), pc_genes)
	for vg_pair in [*variant_gene_pairs]:
		variant_id, gene_id = vg_pair.split(':')

		t.write(variant_id + '\t' + gene_id + '\n')

t.close()

# Make variant vcf file
vcf_file = variant_output_stem + 'all_variant.vcf'
make_variant_vcf_file(variant_gene_pair_file, vcf_file)
