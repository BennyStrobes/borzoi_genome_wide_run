import numpy as np
import os
import sys
import pdb
import gzip
import h5py
import time








######################
# Command line args
######################
file_stem = sys.argv[1]


fast_borzoi_h5_file = file_stem + '_results_fast'
t1 = time.time()
with h5py.File(fast_borzoi_h5_file, 'r') as f:
	snp_chrom = f['snp_chrom'][:]
	snp_pos = f['snp_pos'][:]
	snp = np.asarray(f['snp'][:]).astype(str)
	gene = np.asarray(f['gene'][:]).astype(str)
	logRef = f['logRef'][:]
	logAlt = f['logAlt'][:]

f.close()
t2 = time.time()

fast_borzoi_dicti = {}
for ii, snp_id in enumerate(snp):
	gene_id = gene[ii]
	logRef_val = logRef[ii]
	logAlt_val = logAlt[ii]

	log_sed = logAlt_val - logRef_val
	vg_pair = snp_id + ':' + gene_id
	if vg_pair in fast_borzoi_dicti:
		print('assumption eroror')
		pdb.set_trace()
	fast_borzoi_dicti[vg_pair] = log_sed


'''
# First load in fast borzoi results
fast_borzoi_results_file = file_stem + '_results.txt.gz'
f = gzip.open(fast_borzoi_results_file,'rt')
fast_borzoi_dicti = {}
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	snp_id = data[2]
	gene_id = data[3]
	log_ref = np.asarray(data[4].split(';')).astype(float)
	log_alt = np.asarray(data[5].split(';')).astype(float)
	log_sed = log_alt - log_ref
	vg_pair = snp_id + ':' + gene_id
	if vg_pair in fast_borzoi_dicti:
		print('assumption eroror')
		pdb.set_trace()
	fast_borzoi_dicti[vg_pair] = log_sed
f.close()
'''

sed_h5 = file_stem + '_/sed.h5'
ff = h5py.File(sed_h5, "r")


genes = np.asarray(ff['gene']).astype(str)
snp_indices = np.asarray(ff['si']).astype(int)
snp_ids = (np.asarray(ff['snp'])[snp_indices]).astype(str)
log_seds = np.asarray(ff['logSED']).astype(float)

aa = []
bb = []
for ii, gene_id in enumerate(genes):
	snp_id = snp_ids[ii]
	snp_gene_pair = snp_id + ':' + gene_id

	if snp_gene_pair not in fast_borzoi_dicti:
		print('assumption oeroor')
		pdb.set_trace()

	fast_log_seds = fast_borzoi_dicti[snp_gene_pair]
	norm_log_seds = log_seds[ii]
	aa.append(fast_log_seds)
	bb.append(norm_log_seds)
ff.close()


fast = np.asarray(aa).flatten()
orig = np.asarray(bb).flatten()

print(np.corrcoef(fast,orig)[0,1])

pdb.set_trace()


