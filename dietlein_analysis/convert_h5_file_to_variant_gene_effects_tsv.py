import sys

import h5py
import numpy as np
import pandas as pd

ROW_CHUNK_SIZE = 1024


def load_target_names_and_h5_columns(borzoi_target_file, n_h5_targets):
	targets_df = pd.read_csv(borzoi_target_file, sep="\t", index_col=0)
	identifier_column = None
	for column_name in ("identifier", "target_identifier"):
		if column_name in targets_df.columns:
			identifier_column = column_name
	if identifier_column is not None:
		names = targets_df[identifier_column].astype(str).to_numpy()
	else:
		names = targets_df.iloc[:, 0].astype(str).to_numpy()
	description_column = None
	for column_name in ("description", "target_description"):
		if column_name in targets_df.columns:
			description_column = column_name
	if description_column is not None:
		descriptions = targets_df[description_column].astype(str).to_numpy()
		names = [name + "_" + description for name, description in zip(names, descriptions)]
	names = [name.replace(" ", "_") for name in names]

	if targets_df.shape[0] == n_h5_targets:
		# H5 file already subset to these targets; columns follow target-file row order
		h5_column_indices = np.arange(n_h5_targets)
	else:
		# Target file indexes a subset of the H5 file's columns
		h5_column_indices = np.asarray(targets_df.index, dtype=int)
		if np.any(h5_column_indices < 0) or np.any(h5_column_indices >= n_h5_targets):
			raise ValueError("Target file row count does not match the H5 file and its indices are out of bounds for the H5 file.")

	# De-duplicate repeated target names so tsv columns are unique
	name_counts = {}
	unique_names = []
	for name in names:
		if name in name_counts:
			name_counts[name] = name_counts[name] + 1
			unique_names.append(name + "_" + str(name_counts[name]))
		else:
			name_counts[name] = 0
			unique_names.append(name)

	return unique_names, h5_column_indices


def to_str_array(values):
	return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def main():
	borzoi_h5_file = sys.argv[1]
	borzoi_target_file = sys.argv[2]
	if len(sys.argv) > 3:
		output_tsv_file = sys.argv[3]
	else:
		file_root = borzoi_h5_file[:-3] if borzoi_h5_file.endswith(".h5") else borzoi_h5_file
		output_tsv_file = file_root + "_variant_gene_effects.tsv"

	with h5py.File(borzoi_h5_file, "r") as h5:
		for dataset_name in ("snp", "snp_chrom", "snp_pos", "gene", "logRef", "logAlt"):
			if dataset_name not in h5:
				raise ValueError("Expected dataset '" + dataset_name + "' in the Borzoi H5 file.")

		n_pairs, n_h5_targets = h5["logRef"].shape
		target_names, h5_column_indices = load_target_names_and_h5_columns(borzoi_target_file, n_h5_targets)

		header_columns = ["variant_id", "chrom", "pos", "gene_id"]
		header_columns = header_columns + [name + "_logSED" for name in target_names]
		header_columns = header_columns + [name + "_logRef" for name in target_names]
		header_columns = header_columns + [name + "_logAlt" for name in target_names]

		with open(output_tsv_file, "w") as tsv:
			tsv.write("\t".join(header_columns) + "\n")
			for row_start in range(0, n_pairs, ROW_CHUNK_SIZE):
				row_end = min(row_start + ROW_CHUNK_SIZE, n_pairs)
				variant_ids = to_str_array(h5["snp"][row_start:row_end])
				chroms = to_str_array(h5["snp_chrom"][row_start:row_end])
				positions = h5["snp_pos"][row_start:row_end]
				gene_ids = to_str_array(h5["gene"][row_start:row_end])
				log_ref = h5["logRef"][row_start:row_end, :][:, h5_column_indices]
				log_alt = h5["logAlt"][row_start:row_end, :][:, h5_column_indices]
				log_sed = log_alt - log_ref
				for row_index in range(row_end - row_start):
					fields = [variant_ids[row_index], chroms[row_index], str(int(positions[row_index])), gene_ids[row_index]]
					fields = fields + ["%.6f" % value for value in log_sed[row_index]]
					fields = fields + ["%.6f" % value for value in log_ref[row_index]]
					fields = fields + ["%.6f" % value for value in log_alt[row_index]]
					tsv.write("\t".join(fields) + "\n")

	print("Wrote " + str(n_pairs) + " variant-gene pairs to " + output_tsv_file)


if __name__ == "__main__":
	main()
