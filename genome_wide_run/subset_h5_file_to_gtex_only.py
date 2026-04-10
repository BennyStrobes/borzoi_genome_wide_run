import sys

import h5py
import numpy as np

TARGET_VECTOR_DATASETS = {"logRef", "logAlt"}
ROW_CHUNK_SIZE = 1024


def load_gtex_target_indices(borzoi_gtex_only_target_file):
	return np.loadtxt(
		borzoi_gtex_only_target_file,
		dtype=int,
		delimiter="\t",
		skiprows=1,
		usecols=0,
	)


def dataset_creation_kwargs(dataset):
	kwargs = {"dtype": dataset.dtype}
	if dataset.ndim > 0:
		kwargs["compression"] = dataset.compression
		if dataset.compression is not None:
			kwargs["compression_opts"] = dataset.compression_opts
		if dataset.shuffle:
			kwargs["shuffle"] = dataset.shuffle
		if dataset.fletcher32:
			kwargs["fletcher32"] = dataset.fletcher32
		if dataset.chunks is not None:
			kwargs["chunks"] = dataset.chunks
		if dataset.maxshape is not None:
			kwargs["maxshape"] = dataset.maxshape
	return kwargs


def copy_dataset(source_dataset, destination_h5, gtex_target_indices):
	dataset_name = source_dataset.name.split("/")[-1]
	kwargs = dataset_creation_kwargs(source_dataset)
	if dataset_name in TARGET_VECTOR_DATASETS:
		output_shape = (source_dataset.shape[0], len(gtex_target_indices))
		if "maxshape" in kwargs and kwargs["maxshape"] is not None:
			maxshape = list(kwargs["maxshape"])
			if len(maxshape) > 1:
				maxshape[1] = len(gtex_target_indices)
			kwargs["maxshape"] = tuple(maxshape)
		if "chunks" in kwargs and kwargs["chunks"] is not None:
			chunk_shape = list(kwargs["chunks"])
			chunk_shape[0] = min(chunk_shape[0], max(output_shape[0], 1))
			chunk_shape[1] = min(chunk_shape[1], max(output_shape[1], 1))
			kwargs["chunks"] = tuple(chunk_shape)
		destination_dataset = destination_h5.create_dataset(
			dataset_name,
			shape=output_shape,
			**kwargs,
		)
		for row_start in range(0, source_dataset.shape[0], ROW_CHUNK_SIZE):
			row_end = min(row_start + ROW_CHUNK_SIZE, source_dataset.shape[0])
			destination_dataset[row_start:row_end, :] = source_dataset[row_start:row_end, gtex_target_indices]
	else:
		data = source_dataset[()]
		if "chunks" in kwargs and data.ndim > 0 and kwargs["chunks"] is not None:
			chunk_shape = list(kwargs["chunks"])
			for dim_index, dim_size in enumerate(data.shape):
				chunk_shape[dim_index] = min(chunk_shape[dim_index], dim_size) if dim_size > 0 else 1
			kwargs["chunks"] = tuple(chunk_shape)
		destination_dataset = destination_h5.create_dataset(
			dataset_name,
			data=data,
			**kwargs,
		)
	for attr_key, attr_value in source_dataset.attrs.items():
		destination_dataset.attrs[attr_key] = attr_value


def main():
	borzoi_full_h5_file = sys.argv[1]
	borzoi_gtex_only_h5_file = sys.argv[2]
	borzoi_gtex_only_target_file = sys.argv[3]

	gtex_target_indices = np.atleast_1d(load_gtex_target_indices(borzoi_gtex_only_target_file))

	with h5py.File(borzoi_full_h5_file, "r") as source_h5:
		if "logRef" not in source_h5 or "logAlt" not in source_h5:
			raise ValueError("Expected datasets 'logRef' and 'logAlt' in the full Borzoi H5 file.")

		n_targets = source_h5["logRef"].shape[1]
		if np.any(gtex_target_indices < 0) or np.any(gtex_target_indices >= n_targets):
			raise ValueError("One or more GTEx target indices are out of bounds for the full Borzoi H5 file.")


		with h5py.File(borzoi_gtex_only_h5_file, "w") as destination_h5:
			for attr_key, attr_value in source_h5.attrs.items():
				destination_h5.attrs[attr_key] = attr_value

			for dataset_name in source_h5.keys():
				copy_dataset(source_h5[dataset_name], destination_h5, gtex_target_indices)


if __name__ == "__main__":
	main()
