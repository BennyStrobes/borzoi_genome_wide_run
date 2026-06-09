import numpy as np
import os
import sys
import pdb
import re

def to_underscore(s: str) -> str:
	s = s.replace(")", "")
	# replace " (" OR any run of spaces and hyphens with "_"
	return re.sub(r"(?:\s*\(|[ -]+)", "_", s)

def extract_gtex_tissue_names(gtex_tissue_names_file):
	f = open(gtex_tissue_names_file)
	arr = []
	dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(line)
		dicti[line] = 1

	f.close()
	return np.asarray(arr), dicti

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


def create_mapping_from_gtex_sample_id_to_individual_tissue_format(gtex_sample_attributes_file):
	gtex_sample_to_tissue_name = {}
	gtex_tissues = {}
	f = open(gtex_sample_attributes_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gtex_sample_id = data[0]
		tissue_type = to_underscore(data[6])
		if tissue_type == 'Cells_EBV_transformed_lymphocytes':
			tissue_type = 'Cells_EBV-transformed_lymphocytes'
		elif tissue_type == 'Brain_Spinal_cord_cervical_c_1':
			tissue_type = 'Brain_Spinal_cord_cervical_c-1'


		individual_id = gtex_sample_id.split('-')[0] + '-' + gtex_sample_id.split('-')[1]
		gtex_sample_to_tissue_name[gtex_sample_id] = individual_id + ':' + tissue_type
		gtex_tissues[tissue_type] = 1
	f.close()




	return gtex_sample_to_tissue_name

def create_mapping_from_gtex_tissue_to_target_indices(target_summary_file):
	f = open(target_summary_file)
	dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count= head_count + 1
			continue
		if data[-1] == 'NA':
			continue
		tissue_name = data[-1]
		if tissue_name not in dicti:
			dicti[tissue_name] = []
		dicti[tissue_name].append(int(data[0]))
	f.close()
	return dicti

######################
# Command line args
######################
borzoi_target_file = sys.argv[1]
borzoi_non_gtex_target_file = sys.argv[2]



# Extract borzoi fields
borzoi_target_identifiers, borzoi_target_descriptions = extract_borzoi_target_names(borzoi_target_file)



relevant_borzoi_targets = {
	# Bone / cartilage / musculoskeletal development
	"ENCFF953CRN": "osteoblast",
	"ENCFF092QIW": "osteocyte",
	"ENCFF314QIG": "knee_articular_chondrocyte",
	"ENCFF958WEW": "chondrocyte",
	"ENCFF285OHV": "fetal_skeletal_muscle",
	"ENCFF804YRV": "skeletal_muscle_satellite_cell",

	# Airway / lung / pulmonary development
	"ENCFF153YEN": "bronchial_epithelial_cell",
	"ENCFF383NTW": "bronchial_smooth_muscle_cell",
	"ENCFF734VZN": "alveolar_epithelial_cell",
	"ENCFF595CZA": "airway_epithelial_cell",
	"ENCFF383RPD": "bronchus_fibroblast",
	"ENCFF892OBT": "fetal_lung",

	# Endocrine pancreas / metabolic
	"ENCFF995AUL": "pancreatic_beta_cell",
	"ENCFF225CYB": "endocrine_pancreas_progenitor_cell",

	# Vascular / BP / cardiovascular
	"ENCFF281BWX": "aortic_smooth_muscle_cell",
	"ENCFF537AIY": "coronary_artery_smooth_muscle_cell",
	"ENCFF226UWU": "coronary_artery_endothelial_cell",
	"ENCFF367PUX": "fetal_metanephros",

	# Immune / inflammatory skin
	"ENCFF772RIQ": "regulatory_t_cell",
	"ENCFF362PNM": "activated_cd4_t_cell",
	"ENCFF703WGS": "keratinocyte",

	# Neurodevelopment / neuroendocrine / behavior
	"ENCFF355IUO": "fetal_diencephalon",
	"ENCFF217HQN": "fetal_frontal_cortex",
	"ENCFF674RUW": "ganglionic_eminence_neurosphere",

	# Pigment / hair
	"ENCFF230ZSG": "adult_skin_melanocyte",
	"ENCFF798NEI": "hair_follicle_dermal_papilla_cell",

	# Kidney cell-type / developmental contexts
	"ENCFF863JIL": "glomerular_endothelial_cell",
	"ENCFF346QDJ": "mesangial_cell",
	"ENCFF767MLU": "proximal_tubule_epithelial_cell",

	# Placental / reproductive-development axis
	"ENCFF591NKN": "placental_epithelial_cell",
	"ENCFF779FMR": "placental_pericyte",
	"ENCFF649CUP": "villous_mesenchyme_fibroblast",

	# Additional neurodevelopment / progenitor contexts
	"ENCFF403DKN": "neural_progenitor_cell",
	"ENCFF615IJV": "cortex_neurosphere",
	"ENCFF024MEB": "neural_crest_cell",

	# Skin / epithelial barrier / hair follicle
	"ENCFF343FWG": "dermal_fibroblast",
	"ENCFF544MJM": "hair_follicular_keratinocyte",

	# Vascular / cardiac fibroblast contexts
	"ENCFF256ZZS": "aortic_adventitial_fibroblast",
	"ENCFF287LRZ": "cardiac_ventricle_fibroblast",

	# Additional immune activation states
	"ENCFF122XJD": "activated_t_cell_il2_cd3_cd28",
	"ENCFF330EUN": "activated_cd4_memory_t_cell",
}




t = open(borzoi_non_gtex_target_file,'w')
t.write('orig_target_index\trna_only_target_index\ttarget_identifier\ttarget_description\tGTEx_tissue\n')
new_index = 0
for ii, identifier in enumerate(borzoi_target_identifiers):

	if identifier.startswith('GTEX'):
		continue

	if borzoi_target_descriptions[ii].startswith('RNA') == False:
		continue
	if identifier not in relevant_borzoi_targets:
		continue
	t.write(str(ii) + '\t' + str(new_index) + '\t' + identifier + '\t' + borzoi_target_descriptions[ii] + '\t' + relevant_borzoi_targets[identifier] + '\n')
	new_index = new_index + 1
t.close()




