import sys


def load_specified_gene_ids(comma_separated_gene_ids):
	# Match on un-versioned ensembl ids so version suffixes don't have to agree
	gene_ids = {}
	for gene_id in comma_separated_gene_ids.split(","):
		gene_ids[gene_id.split(".")[0]] = 1
	return gene_ids


def main():
	input_tsv_file = sys.argv[1]
	output_tsv_file = sys.argv[2]
	specified_gene_ids = load_specified_gene_ids(sys.argv[3])

	n_kept = 0
	f = open(input_tsv_file)
	t = open(output_tsv_file, "w")
	head_count = 0
	for line in f:
		line = line.rstrip("\n")
		data = line.split("\t")
		if head_count == 0:
			head_count = 1
			if "gene_id" not in data:
				raise ValueError("Expected 'gene_id' column in " + input_tsv_file)
			gene_id_column = data.index("gene_id")
			t.write(line + "\n")
			continue
		if data[gene_id_column].split(".")[0] in specified_gene_ids:
			t.write(line + "\n")
			n_kept = n_kept + 1
	f.close()
	t.close()

	print("Wrote " + str(n_kept) + " variant-gene pairs to " + output_tsv_file)


if __name__ == "__main__":
	main()
