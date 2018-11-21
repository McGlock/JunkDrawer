import pandas as pd
import sys
from os import path


# Uses prot.accession2taxid.lineage.tsv created by acc_lineage_mapper.py
#     to build accession_id_lineage_map.tsv file required by TreeSAPP
#     in the building of a refpkg.
# 1. acc_lineage_mapper.py output file with acc and lineage info
# 2. Cleaned fasta file with NCBI accession as header name
# 3. Output path for file used by TreeSAPP
#    *must be named 'accession_id_lineage_map.tsv'
# Usage:
# python build_ref_lineage.py [1] [2] [3]

acc_lin_map_file = sys.argv[1]
clean_fasta_file = sys.argv[2]
output_path = sys.argv[3]
output_file = path.join(output_path, 'accession_id_lineage_map.tsv')
exclude_acc_file = path.join(output_path, clean_fasta_file.split('.')[0] + '.exclude.txt')

acc_lin_map_chunk = pd.read_csv(acc_lin_map_file, sep='\t', header=0, chunksize=50000000)
acc_list = []

with open(clean_fasta_file, 'r') as i:
    data = i.readlines()
    for line in data:
        if '>' in line:
            acc_id = line.split(' ')[0].strip('>')
            acc_list.append(acc_id)

df_list = []
for chunk in acc_lin_map_chunk:
	filter_df = chunk[chunk['accession.version'].isin(acc_list)]
	filter_df['full_lineage'] = [x + y for x, y in
								zip(filter_df['lineage'], filter_df['sciname'])]
	sub_df = filter_df[['accession.version', 'full_lineage']]
	#sub_df['accession.version'] = [x.replace('_', '') for x in sub_df['accession.version']]
	df_list.append(sub_df)
combine_df = pd.concat(df_list)
combine_df.to_csv(output_file, sep='\t', index=False, header=False)

exclude_acc_list = []
for acc in acc_list:
	if acc not in list(combine_df['accession.version']):
		exclude_acc_list.append(acc)

with open(exclude_acc_file, 'w') as e:
	e.write('\n'.join(exclude_acc_list))