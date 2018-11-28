import pandas as pd
import sys


# Maps NCBI accession.version IDs to NCBI Taxonomy IDs and lineages
# Outputs a merged table containing all the info from the following files.
# 1. prot.accession2taxid file downloaded from NCBI Taxonony FTP
# 2. fullnamelineage.dmp file from downloaded NCBI Taxonomy dump
# 3. prot.accession2taxid.lineage.tsv or other output filename
# 4. merged.dmp file from downloaded NCBI Taxonomy dump
# Useage: 
# python acc_lineage_mapper.py [1] [2] [3]

acc2taxid_file = sys.argv[1]
taxid2lin_file = sys.argv[2]
output_file = sys.argv[3]
merged_taxa_file = sys.argv[4]


taxid2lin_df = pd.read_csv(taxid2lin_file, sep='\t\\|\t',
						names=['taxid', 'sciname', 'lineage'], engine='python')
taxid2lin_df['lineage'] = [x.strip('\t|') for x in taxid2lin_df['lineage']]

merged_taxa_df = pd.read_csv(merged_taxa_file, sep='\t\\|\t',
						names=['taxid', 'merge2taxid'], engine='python')
merged_taxa_df['merge2taxid'] = [int(x.strip('\t|')) for x in merged_taxa_df['merge2taxid']]

merged_taxa_dict = {x[0]:x[1] for x in
						zip(merged_taxa_df['taxid'], merged_taxa_df['merge2taxid'])
						}

with open(output_file, 'w') as f:
	f.write('\t'.join(['accession', 'accession.version', 'taxid', 'gi',
						'sciname', 'lineage\n'])) # This assumes column names for now

acc2taxid_chunks = pd.read_table(acc2taxid_file, sep='\t', header=0, chunksize=50000000)

with open(output_file, 'a') as a:
	for count, chunk in enumerate(acc2taxid_chunks, start=1):	
		chunk['taxid'] = [merged_taxa_dict[x] if x in merged_taxa_dict.keys()
							else x	for x in chunk['taxid']
							]
		comb_df = pd.merge(chunk, taxid2lin_df, on='taxid', how='left')
		comb_df.to_csv(a, sep='\t', header=False, index=False)
		print('Chunk {} mapped to lineage info...'.format(count))