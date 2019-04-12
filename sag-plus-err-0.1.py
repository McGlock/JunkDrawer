import matplotlib
matplotlib.use('agg')
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys
from os.path import join as joinpath
from functools import reduce
import numpy as np
from collections import Counter


def calc_err(df):
	# build error type df for each filter separately
	group_df = df.groupby(['sag_id', 'algorithm', 'level']).sum().reset_index()
	group_df['precision'] = group_df['TruePos'] / \
								(group_df['TruePos'] + group_df['FalsePos'])
	group_df['sensitivity'] = group_df['TruePos'] / \
								(group_df['TruePos'] + group_df['FalseNeg'])
	group_df['specificity'] = group_df['TrueNeg'] / \
								(group_df['TrueNeg'] + group_df['FalsePos'])
	group_df['type1_error'] = group_df['FalsePos'] / \
								(group_df['FalsePos'] + group_df['TrueNeg'])
	group_df['type2_error'] = group_df['FalseNeg'] / \
								(group_df['FalseNeg'] + group_df['TruePos'])
	group_df['F1_score'] = 2 * ((group_df['precision'] * group_df['sensitivity']) / \
								(group_df['precision'] + group_df['sensitivity']))
	group_df.set_index(['sag_id', 'algorithm', 'level'], inplace=True)
	stats_df = group_df[['precision', 'sensitivity', 'specificity', 'type1_error',
							'type2_error', 'F1_score']]
	stack_df = stats_df.stack().reset_index()
	stack_df.columns = ['sag_id', 'algorithm', 'level', 'statistic', 'score']
	return stack_df

# Map genome id and contig id to taxid for error analysis
sag_tax_map = '/home/rmclaughlin/Ryan/CAMI_gold/CAMI_I_HIGH/genome_taxa_info.tsv'
sag_taxmap_df = pd.read_csv(sag_tax_map, sep='\t', header=0)
sag_taxmap_df['sp_taxid'] = [int(x) for x in sag_taxmap_df['@@TAXID']]
sag_taxmap_df['sp_name'] = [x.split('|')[-2] for x in sag_taxmap_df['TAXPATHSN']]
taxpath_list = [[str(x) for x in x.split('.')[0].split('|')]
					for x in sag_taxmap_df['TAXPATH']
					]
taxpath_df = pd.DataFrame(taxpath_list, columns=['domain', 'phylum', 'class', 'order',
													'family', 'genus', 'species', 'strain'
													])
taxpath_df['CAMI_genomeID'] = [x for x in sag_taxmap_df['_CAMI_genomeID']]
# fix empty species id's
taxpath_df['species'] = [x[1] if str(x[0]) == '' else x[0] for x in 
							zip(taxpath_df['species'], taxpath_df['genus'])
							]
# Map MetaG contigs to their genomes
mg_contig_map = '/home/rmclaughlin/Ryan/CAMI_gold/CAMI_I_HIGH/' + \
				'gsa_mapping_pool.binning.trimmed'
mg_contig_map_df = pd.read_csv(mg_contig_map, sep='\t', header=0)
mg_contig_map_df['TAXID'] = [str(x) for x in mg_contig_map_df['TAXID']]

# Merge contig map and taxpath DFs
tax_mg_df = taxpath_df.merge(mg_contig_map_df, left_on='strain', right_on='TAXID',
								how='right'
								)
tax_mg_df = tax_mg_df[['@@SEQUENCEID', 'CAMI_genomeID', 'domain', 'phylum', 'class', 'order',
						'family', 'genus', 'species', 'strain'
						]]
files_path = '/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/'

# Build error df from recruits files
# MinHash
mh_path = joinpath(files_path, 'minhash_recruits/')
mh_df_list = []
mh_file_list = [x for x in os.listdir(mh_path) 
					if 'mhr_recruits.tsv' in x
					]
for mh_file in mh_file_list:
	file_path = os.path.join(mh_path, mh_file)
	file_df = pd.read_csv(file_path, sep='\t', header=None,
							names=['sag_id', 'subcontig_id', 'contig_id']
							)
	mh_df_list.append(file_df)
mh_concat_df = pd.concat(mh_df_list)

# RPKM
rpkm_path = joinpath(files_path, 'rpkm_recruits/')
rpkm_df_list = []
rpkm_file_list = [x for x in os.listdir(rpkm_path) 
					if 'ara_recruits.tsv' in x
					]
for rpkm_file in rpkm_file_list:
	file_path = os.path.join(rpkm_path, rpkm_file)
	file_df = pd.read_csv(file_path, sep='\t', header=None,
							names=['sag_id', 'contig_id']
							)
	rpkm_df_list.append(file_df)
rpkm_concat_df = pd.concat(rpkm_df_list)
rpkm_concat_df.insert(1, 'subcontig_id', '')

# Tetra GMM
tetra_path = joinpath(files_path, 'tetra_recruits/')
tetra_df_list = []
tetra_file_list = [x for x in os.listdir(tetra_path) 
					if 'tra_recruits.tsv' in x
					]
for tetra_file in tetra_file_list:
	file_path = os.path.join(tetra_path, tetra_file)
	file_df = pd.read_csv(file_path, sep='\t', header=None,
							names=['sag_id', 'subcontig_id', 'contig_id']
							)
	tetra_df_list.append(file_df)
tetra_concat_df = pd.concat(tetra_df_list)

# Final Recruits
final_file = joinpath(files_path, 'final_recruits/final_recruits.tsv')
final_df = pd.read_csv(final_file, sep='\t', header=0, index_col=0,
							names=['sag_id', 'subcontig_id', 'contig_id']
							)

mh_concat_df['algorithm'] = 'MinHash'
rpkm_concat_df['algorithm'] = 'RPKM'
tetra_concat_df['algorithm'] = 'tetra_GMM'
final_df['algorithm'] = 'combined'
final_concat_df = pd.concat([mh_concat_df, rpkm_concat_df, tetra_concat_df, final_df])
final_tax_df = final_concat_df.merge(tax_mg_df, left_on='contig_id', right_on='@@SEQUENCEID',
								how='left'
								)

error_list = []
for i, sag_id in enumerate(list(set(final_tax_df['sag_id']))):
	sag_key_list = [str(s) for s in set(tax_mg_df['CAMI_genomeID']) if str(s) in sag_id]
	sag_key = max(sag_key_list, key=len)
	sag_sub_df = final_tax_df.loc[final_tax_df['sag_id'] == sag_id]
	for algo in set(sag_sub_df['algorithm']):
		algo_sub_df = sag_sub_df.loc[sag_sub_df['algorithm'] == algo]
		for col in taxpath_df.columns:
			col_key = final_tax_df.loc[final_tax_df['CAMI_genomeID'] == sag_key,
										col].iloc[0]
			contig2col_key = {x[0]: x[1] for x in zip(tax_mg_df['@@SEQUENCEID'],
												tax_mg_df[col])
												}
			print(i, sag_id, algo, col)
			if col == 'CAMI_genomeID':
				col = 'perfect'
				col_key = sag_key
			for contig_id in contig2col_key.keys():
				err_list = [sag_id, algo, col, contig_id, 0, 0, 0, 0]
				if ((contig_id in list(algo_sub_df['contig_id'])) and 
					(contig2col_key[contig_id] == col_key)
					):
					err_list[4] = 1 # 'TruePos'
				elif ((contig_id in list(algo_sub_df['contig_id'])) and 
					(contig2col_key[contig_id] != col_key)
					):
					err_list[5] = 1 # 'FalsePos'
				elif  ((contig_id not in list(algo_sub_df['contig_id'])) and
					(contig2col_key[contig_id] == col_key)
					):
					err_list[6] = 1 # 'FalseNeg'
				elif ((contig_id not in list(algo_sub_df['contig_id'])) and 
					(contig2col_key[contig_id] != col_key)
					):
					err_list[7] = 1 # 'TrueNeg'
				error_list.append(err_list)

final_err_df = pd.DataFrame(error_list, columns=['sag_id', 'algorithm', 'level',
													'contig_id', 'TruePos', 'FalsePos',
													'FalseNeg', 'TrueNeg'
													])

#final_err_df.to_csv('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/' + \
#					'error_analysis/All_error.tsv', index=False, sep='\t'
#					)

calc_stats_df = calc_err(final_err_df)
calc_stats_df.to_csv('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/' + \
					'error_analysis/All_stats.tsv', index=False, sep='\t'
					)
for level in set(calc_stats_df['level']):
	level_df = calc_stats_df.loc[calc_stats_df['level'] == level]
	sns.set_context("paper")
	ax = sns.catplot(x="statistic", y="score", hue='algorithm', kind='box',
						data=level_df, aspect=2
						)
	plt.title('SAG-plus CAMI-1-High error analysis')
	ax._legend.set_title('Filter Type')
	plt.savefig('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/' + \
				'error_analysis/' + level + '_error_boxplox.svg',
				bbox_inches='tight'
				)
	plt.clf()


'''



# Recruited fastas
fasta_path = joinpath(files_path, 'final_recruits/')
fasta_file_list = []
fasta_file_list = [x for x in os.listdir(fasta_path) 
					if 'final_recruits.fasta' in x
					]
sag_taxid_list = []
for fa_file in fasta_file_list:
	sag_id = fa_file.split('.final_recruits.fasta')[0]
	print(sag_id)
	file_path = os.path.join(fasta_path, fa_file)
	with open(file_path, 'r') as fast_in:
		data = fast_in.readlines()
		for line in data:
			if '>' in line:
				sag_key_list = [s for s in sag_taxmap_df['_CAMI_genomeID'] if s in sag_id]
				sag_key = max(sag_key_list, key=len)
				taxid = str(sag_taxmap_df.loc[sag_taxmap_df['_CAMI_genomeID'] == sag_key,
											'@@TAXID'].values[0]).split('.')[0]
				mg_contid_id = line.rsplit('_',1)[0].lstrip('>')
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
		cami_include_ids = list(taxpath_df.loc[taxpath_df['CAMI_genomeID'] == sag_key,
								'CAMI_genomeID']
								)
		mg_tot_contigs = mg_contig_map_df.loc[mg_contig_map_df['BINID'
									].isin(cami_include_ids)]['@@SEQUENCEID'].count()
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
sag_precision_df.to_csv('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/' + \
						'error_analysis/stats_level.tsv', index=False, sep='\t'
						)
print(len(set(sag_precision_df['sag_id'])))

# build multi-level precision boxplot
sns.set_context("paper")
ax = sns.catplot( x="level", y="score", hue='statistic', kind='box',
					data=sag_precision_df, aspect=2
					)
plt.title('SAG-plus CAMI-1-High')
plt.savefig('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/' + \
			'error_analysis/multi-level_precision_boxplox.svg', bbox_inches='tight'
			)
plt.clf()
'''