import numpy as np
import os
import sys
import pdb






#######################
# Command line args
#######################
fm_summary_file = sys.argv[1]
fm_variant_vcf_file = sys.argv[2]



# Extract dictionary list of variants
dicti = {}
f = open(fm_summary_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	dicti[data[0]] = 1
f.close()



t = open(fm_variant_vcf_file,'w')
for variant_id in np.sort([*dicti]):
	var_info = variant_id.split('_')
	t.write(var_info[0] + '\t' + var_info[1] + '\t' + variant_id + '\t' + var_info[2] + '\t' + var_info[3] + '\n')
t.close()