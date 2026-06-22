import numpy as np
import os
import sys
import pdb
import re


def extract_borzoi_target_names(borzoi_target_file):
	f = open(borzoi_target_file)
	head_count = 0
	target_identifiers = []
	target_descriptions = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		identifier = data[1]
		if identifier.endswith('+'):
			identifier = identifier.split('+')[0]
		elif identifier.endswith('-'):
			continue
		description = data[-1]
		target_identifiers.append(identifier)
		target_descriptions.append(description)
	f.close()

	return np.asarray(target_identifiers), np.asarray(target_descriptions)


######################
# Command line args
######################
borzoi_target_file = sys.argv[1]
borzoi_all_target_file = sys.argv[2]


# Extract borzoi fields (+ strand only; '-' strand rows are skipped so that the
# row index lines up with the target axis of the prediction H5 files)
borzoi_target_identifiers, borzoi_target_descriptions = extract_borzoi_target_names(borzoi_target_file)


t = open(borzoi_all_target_file, 'w')
t.write('orig_target_index\tall_target_index\ttarget_identifier\ttarget_description\tassay\n')
new_index = 0
for ii, identifier in enumerate(borzoi_target_identifiers):
	description = borzoi_target_descriptions[ii]
	assay = description.split(':')[0] if ':' in description else 'NA'
	t.write(str(ii) + '\t' + str(new_index) + '\t' + identifier + '\t' + description + '\t' + assay + '\n')
	new_index = new_index + 1
t.close()
