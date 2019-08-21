import matplotlib
matplotlib.use('agg')
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys


# Map genome id and contig id to taxid for error analysis
sag_tax_map = '/home/rmclaughlin/Ryan/CAMI_gold/CAMI_I_HIGH/genome_taxa_info.tsv'
sag_taxmap_df = pd.read_csv(sag_tax_map, sep='\t', header=0)
sag_taxmap_df['sp_taxid'] = [int(x) for x in sag_taxmap_df['@@TAXID']]
sag_taxmap_df['sp_name'] = [x.split('|')[-2] for x in sag_taxmap_df['TAXPATHSN']]
taxpath_list = [x.split('.')[0].split('|') for x in sag_taxmap_df['TAXPATH']]
taxpath_df = pd.DataFrame(taxpath_list, columns=['domain', 'phylum', 'class', 'order',
													'family', 'genus', 'species', 'strain'
													])
taxpath_df['CAMI_genomeID'] = [x for x in sag_taxmap_df['_CAMI_genomeID']]
# fix empty species id's
taxpath_df['species'] = [x[1] if str(x[0]) == '' else x[0] for x in 
							zip(taxpath_df['species'], taxpath_df['genus'])
							]

# Map MetaG contigs to their genomes
mg_contig_map = '/home/rmclaughlin/Ryan/CAMI_gold/CAMI_I_HIGH/gsa_mapping_pool.binning.trimmed'
mg_contig_map_df = pd.read_csv(mg_contig_map, sep='\t', header=0)

files_path = '/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/test_sourmash/'
df_list = []
sag_taxid_list = []
# only process if error file exists
error_file_list = [x for x in os.listdir(files_path) 
					if ('_error_stats.tsv' in x and
					x != 'total_error_stats.tsv')
					]
fasta_file_list = [x.rsplit('_', 2)[0] + '.idf_plus_htf.predicted_contigs.fasta'
					for x in error_file_list
					]

for err_file in error_file_list:
	file_path = os.path.join(files_path, err_file)
	taxid_list = []
	file_df = pd.read_csv(file_path, sep='\t', header=0)
	df_list.append(file_df)
	

for fa_file in fasta_file_list:
	sag_id = fa_file.split('.idf_plus_htf.predicted_contigs.fasta')[0]
	file_path = os.path.join(files_path, fa_file)
	with open(file_path, 'r') as fast_in:
		data = fast_in.readlines()
		for line in data:
			if '>' in line:
				taxid = line.split('|')[-1].rstrip('\n')
				mg_contid_id = line.split('_')[0].lstrip('>')
				sag_key_list = [s for s in sag_taxmap_df['_CAMI_genomeID'] if s in sag_id]
				sag_key = max(sag_key_list, key=len)
				sag_taxid_list.append([mg_contid_id, sag_id, taxid])

sag_taxid_df = pd.DataFrame(sag_taxid_list, columns=['mg_contig_id', 'sag_id', 'strain'])
sag_lineage_df = sag_taxid_df.merge(taxpath_df, on='strain', how='left')
sag_lineage_df.drop_duplicates(inplace=True)

# Calculate all level of precision
sag_precision_list = []
for sag_id in list(set(sag_lineage_df['sag_id'])):
	sag_key_list = [s for s in sag_taxmap_df['_CAMI_genomeID'] if s in sag_id]
	sag_key = max(sag_key_list, key=len)
	sag_sub_df = sag_lineage_df.loc[sag_lineage_df['sag_id'] == sag_id]

	for col in taxpath_df.columns:
		col_key = taxpath_df.loc[taxpath_df['CAMI_genomeID'] == sag_key, col].iloc[0]
		cami_include_ids = list(taxpath_df.loc[taxpath_df['CAMI_genomeID'] == sag_key, 'CAMI_genomeID'])
		mg_tot_contigs = mg_contig_map_df.loc[mg_contig_map_df['BINID'].isin(cami_include_ids)]['@@SEQUENCEID'].count()

		match_dict = sag_sub_df[col].value_counts().to_dict()
		if col == 'CAMI_genomeID':
			col = 'perfect'
		sag_count = int(match_dict[col_key])
		all_count = sum(match_dict.values())
		precision = sag_count/all_count
		sensitivity = sag_count/mg_tot_contigs
		if sensitivity > 1.0:
			sensitivity = 1.0
		F1_score = 2 * ((precision * sensitivity) / (precision + sensitivity))
		#sag_precision_list.append([sag_id, col, precision])
		sag_precision_list.append([sag_id, col, 'precision', precision])
		sag_precision_list.append([sag_id, col, 'sensitivity', sensitivity])
		sag_precision_list.append([sag_id, col, 'F1_score', F1_score])

sag_precision_df = pd.DataFrame(sag_precision_list,
								columns=['sag_id', 'level', 'statistic', 'score']
								)
sag_precision_df.to_csv('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/SAG_plus_error/' \
						'stats_level.tsv', index=False, sep='\t'
						)
print(len(set(sag_precision_df['sag_id'])))

# build multi-level precision boxplot
sns.set_context("paper")
ax = sns.catplot( x="level", y="score", hue='statistic', kind='box', data=sag_precision_df, aspect=2)
plt.title('SAG-plus CAMI-1-High')
plt.savefig('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/SAG_plus_error/multi-level_precision_boxplox.svg', bbox_inches='tight')
plt.clf()


# concat all error dfs
final_err_df = pd.concat(df_list)
final_err_df.to_csv('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/SAG_plus_error/All_error.tsv', index=False, sep='\t')
filter_out = ['TruePos', 'TrueNeg', 'FalsePos', 'FalseNeg']
filter_final_err_df = final_err_df.loc[~final_err_df.filter_type.isin(filter_out)]
# build boxplot for all errors analysis
sns.set_context("paper")
ax = sns.catplot( x="statistic", y="score", hue='filter_type', kind='box', data=filter_final_err_df, aspect=2)
plt.title('SAG-plus CAMI-1-High error analysis')
ax._legend.set_title('Filter Type')
plt.savefig('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/SAG_plus_error/perfect_match_error_boxplox.svg', bbox_inches='tight')
plt.clf()


