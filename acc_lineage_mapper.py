import pandas as pd
import sys


# Maps NCBI accession.version IDs to NCBI Taxonomy IDs and lineages
# Outputs a merged table containing all the info from the following files.
# 1. prot.accession2taxid file downloaded from NCBI Taxonony FTP
# 2. fullnamelineage.dmp file from downloaded NCBI Taxonomy dump
# 3. prot.accession2taxid.lineage.tsv or other output filename
# Useage: 
# python acc_lineage_mapper.py [1] [2] [3]

acc2taxid_file = sys.argv[1]
taxid2lin_file = sys.argv[2]
output_file = sys.argv[3]

acc2taxid_chunks = pd.read_table(acc2taxid_file, sep='\t', header=0, chunksize=50000000)

taxid2lin_df = pd.read_csv(taxid2lin_file, sep='\t\\|\t',
						names=['taxid', 'sciname', 'lineage'], engine='python')
taxid2lin_df['lineage'] = [x.strip('\t|') for x in taxid2lin_df['lineage']]

with open(output_file, 'w') as f:
	f.write('\t'.join(['accession', 'accession.version', 'taxid', 'gi',
						'sciname', 'lineage\n'])) # This assumes column names for now

with open(output_file, 'a') as a:
	for count, chunk in enumerate(acc2taxid_chunks, start=1):
		merge_df = pd.merge(chunk, taxid2lin_df, on='taxid', how='left')
		merge_df.to_csv(a, sep='\t', header=False, index=False)
		print('Chunk {} mapped to lineage info...'.format(count))