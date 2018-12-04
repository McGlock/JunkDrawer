import pandas as pd
import sys
from os import path
from pathlib import Path
import numpy as np


# Uses prot.accession2taxid.lineage.tsv created by acc_lineage_mapper.py
#     to build accession_id_lineage_map.tsv file required by TreeSAPP
#     in the building of a refpkg.
# 1. acc_lineage_mapper.py output file with acc and lineage info
# 2. Cleaned fasta file with NCBI accession as header name
# 3. Output path for file used by TreeSAPP
#    *must be named 'accession_id_lineage_map.tsv'
# Usage:
# python build_ref_lineage.py [1] [2] [3]


def fasta_slurp(fasta_file):
	# assumes header begins with accession followed by a space
	# assumes no accession duplicates
	fa_dict = {}
	with open(fasta_file, 'r') as i:
		data = i.read().split('>')
		for record in data[1:]:
			header, seq = record.split('\n', 1)
			acc_id = header.split(' ')[0]
			fa_dict[acc_id] = (header, seq)
	return fa_dict


# Set file names
acc_lin_map_file = sys.argv[1]
clean_fasta_file = sys.argv[2]
output_path = sys.argv[3]
output_file = path.join(output_path, 'accession_id_lineage_map.tsv')
tax_id_file = clean_fasta_file.split('.')[0] + '.tax_ids.txt'
nolin_acc_file = clean_fasta_file.split('.')[0] + '.nolin.txt'
taxed_fasta_file = clean_fasta_file.split('.')[0] + '.taxed.fasta'
acc_lin_map_custom_file = clean_fasta_file.split('.')[0] + '.accession_id_lineage_map.tsv'
missing_lin_map_custom_file = clean_fasta_file.split('.')[0] + '.missing_lin_map.tsv'

# Open accession table in chunks due to its size
acc_lin_map_chunk = pd.read_csv(acc_lin_map_file, sep='\t', header=0, chunksize=50000000)

# Build fasta dictionary from recruited sequences
fasta_dict = fasta_slurp(clean_fasta_file)

# Build a dataframe of the accession to lineage map for the recruited seqs
df_list = []
nan_list = []
for chunk in acc_lin_map_chunk:
	filter_df = chunk[chunk['accession.version'].isin(fasta_dict.keys())]
	filter_df['full_lineage'] = [x + y for x, y in
								zip(filter_df['lineage'], filter_df['sciname'])]
	no_nan_df = filter_df[filter_df['full_lineage'].notna()]
	nan_df = filter_df[filter_df['full_lineage'].isna()]
	df_list.append(no_nan_df)
	nan_list.append(nan_df)
combine_df = pd.concat(df_list)
if not Path(acc_lin_map_custom_file).is_file():
	combine_df.to_csv(acc_lin_map_custom_file, sep='\t', header=True, index=False)
missing_lin_df = pd.concat(nan_list)
if not Path(missing_lin_map_custom_file).is_file():
	missing_lin_df.to_csv(missing_lin_map_custom_file, sep='\t', header=True, index=False)
sub_df = combine_df[['accession.version', 'full_lineage']]
sub_df.to_csv(output_file, sep='\t', index=False, header=False)
#for x in zip(sub_df['accession.version'], sub_df['full_lineage']):
#	print(x)
#	print(x[1].split(';')[-1] + ' | ' + x[0])
sub_df['spp | acc'] = [x[1].split(';')[-1] + ' | ' + x[0] for x in
						zip(sub_df['accession.version'], sub_df['full_lineage'])
						]
sub_df['lineage2genus'] = [x.rsplit(';', 1)[0] for x in sub_df['full_lineage']]
sub_df = sub_df.reset_index()
tax_id_df = sub_df[['spp | acc', 'lineage2genus']]
tax_id_df.to_csv(tax_id_file, sep='\t', header=False)

# Add taxid to the begining of the seq headers
with open(taxed_fasta_file , 'w') as t:
	for k, v in fasta_dict.items():
		acc_id = k
		header, seq = v
		if acc_id in list(combine_df['accession.version']):
			taxid = str(combine_df[combine_df['accession.version'] == acc_id].iloc[0]['taxid'])
			lineage = combine_df[combine_df['accession.version'] == acc_id
									].iloc[0]['full_lineage']
			new_id = '|'.join(['>' + taxid, acc_id])
			#new_header = ' '.join([new_id, lineage])
			new_record = '\n'.join([new_id, seq])
			t.write(new_record)

# Check if there are any accessions that didn't get lineage map
nolin_acc_list = []
for acc in fasta_dict.keys():
	if acc not in list(combine_df['accession.version']):
		nolin_acc_list.append(acc)

# Save nolin accession so they can be removed
with open(nolin_acc_file, 'w') as e:
	e.write('\n'.join(nolin_acc_list))