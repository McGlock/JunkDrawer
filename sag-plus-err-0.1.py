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
from os import listdir, makedirs, path


def calc_err(df):
	# build error type df for each filter separately
	group_df = df.copy()
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
sag_tax_map = '/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/genome_taxa_info.tsv'
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
mg_contig_map = '/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/' + \
				'gsa_mapping_pool.binning.trimmed'
mg_contig_map_df = pd.read_csv(mg_contig_map, sep='\t', header=0)
mg_contig_map_df['TAXID'] = [str(x) for x in mg_contig_map_df['TAXID']]

# Merge contig map and taxpath DFs
tax_mg_df = taxpath_df.merge(mg_contig_map_df, left_on='CAMI_genomeID', right_on='BINID',
								how='right'
								)
tax_mg_df = tax_mg_df[['@@SEQUENCEID', 'CAMI_genomeID', 'domain', 'phylum', 'class', 'order',
						'family', 'genus', 'species', 'strain'
						]]
files_path = sys.argv[1]
err_path = files_path + '/error_analysis'
if not path.exists(err_path):
	makedirs(err_path)

# Build error df from recruits files
# MinHash
mh_path = joinpath(files_path, 'minhash_recruits/')
mh_df_list = []
mh_file_list = [x for x in os.listdir(mh_path) 
					if 'mhr_recruits.tsv' in x
					]
print('loading minhash files')
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
print('loading rpkm files')
for rpkm_file in rpkm_file_list:
	file_path = os.path.join(rpkm_path, rpkm_file)
	file_df = pd.read_csv(file_path, sep='\t', header=None,
							names=['sag_id', 'subcontig_id', 'contig_id']
							)
	rpkm_df_list.append(file_df)
rpkm_concat_df = pd.concat(rpkm_df_list)

# Tetra GMM
tetra_path = joinpath(files_path, 'tetra_recruits/')
tetra_df_list = []
tetra_file_list = [x for x in os.listdir(tetra_path) 
					if 'tra_recruits.tsv' in x
					]
print('loading tetra files')
for tetra_file in tetra_file_list:
	file_path = os.path.join(tetra_path, tetra_file)
	file_df = pd.read_csv(file_path, sep='\t', header=None,
							names=['sag_id', 'subcontig_id', 'contig_id']
							)
	tetra_df_list.append(file_df)
tetra_concat_df = pd.concat(tetra_df_list)

# Final Recruits
final_file = joinpath(files_path, 'final_recruits/final_recruits.tsv')
print('loading combined files')
final_df = pd.read_csv(final_file, sep='\t', header=0, index_col=0,
							names=['sag_id', 'contig_id']
							)
final_df['subcontig_id'] = None

mh_concat_df['algorithm'] = 'MinHash'
rpkm_concat_df['algorithm'] = 'RPKM'
#final_concat_df = pd.concat([mh_concat_df, rpkm_concat_df])


tetra_concat_df['algorithm'] = 'tetra_GMM'
final_df['algorithm'] = 'combined'
final_concat_df = pd.concat([mh_concat_df, rpkm_concat_df, tetra_concat_df, final_df])
final_group_df = final_concat_df.groupby(['sag_id', 'algorithm', 'contig_id'])[
											'subcontig_id'].count().reset_index()
print(mh_concat_df.head())
print(rpkm_concat_df.head())
print(tetra_concat_df.head())
print(final_df.head())

print('merging all')
final_tax_df = final_group_df.merge(tax_mg_df, left_on='contig_id', right_on='@@SEQUENCEID',
								how='left'
								)
sag_cnt_dict = final_tax_df.groupby('sag_id')['sag_id'].count().to_dict()
#final_tax_df['contig_count'] = [sag_cnt_dict[x] for x in final_tax_df['sag_id']]
#final_tax_df.sort_values('contig_count', ascending=True, inplace=True)
error_list = []
algo_list = ['MinHash', 'RPKM', 'tetra_GMM', 'combined']
level_list = ['genus', 'species', 'strain', 'CAMI_genomeID']
for i, sag_id in enumerate(list(final_df['sag_id'].unique())):
	sag_key_list = [str(s) for s in set(tax_mg_df['CAMI_genomeID']) if str(s) in sag_id]
	sag_key = max(sag_key_list, key=len)
	sag_sub_df = final_tax_df.loc[final_tax_df['sag_id'] == sag_id]
	for algo in algo_list:
		algo_sub_df = sag_sub_df.loc[sag_sub_df['algorithm'] == algo]
		for col in level_list:
			col_key = final_tax_df.loc[final_tax_df['CAMI_genomeID'] == sag_key,
										col].iloc[0]
			cami_include_ids = list(set(tax_mg_df.loc[tax_mg_df[col] == col_key,
									'CAMI_genomeID'])
									)
			mg_include_contigs = list(set(tax_mg_df.loc[tax_mg_df['CAMI_genomeID'
									].isin(cami_include_ids)]['@@SEQUENCEID'])
									)
			sag_include_contigs = list(set(tax_mg_df.loc[tax_mg_df['CAMI_genomeID'
									].isin([sag_key])]['@@SEQUENCEID'])
									)

			print(i, sag_id, algo, col, col_key, len(mg_include_contigs),
				len(sag_include_contigs), len(cami_include_ids)
				)
			if col == 'CAMI_genomeID':
				col = 'perfect'
				col_key = sag_key
			err_list = [sag_id, algo, col, 0, 0, 0, 0]
			for contig_id in tax_mg_df['@@SEQUENCEID']:
				if contig_id in list(algo_sub_df['contig_id']):
					if contig_id in mg_include_contigs:
						err_list[3] += 1 # 'TruePos'
					else:
						err_list[4] += 1 # 'FalsePos'
				else:
					if contig_id in sag_include_contigs:
						err_list[5] += 1 # 'FalseNeg'
					else:
						err_list[6] += 1 # 'TrueNeg'
			error_list.append(err_list)

final_err_df = pd.DataFrame(error_list, columns=['sag_id', 'algorithm', 'level',
													'TruePos', 'FalsePos',
													'FalseNeg', 'TrueNeg'
													])

final_err_df.to_csv(err_path + '/All_error.tsv', index=False, sep='\t')

calc_stats_df = calc_err(final_err_df)
stat_list = ['precision', 'sensitivity', 'F1_score']
calc_stats_df = calc_stats_df.loc[calc_stats_df['statistic'].isin(stat_list)]
calc_stats_df.to_csv(err_path + '/All_stats.tsv', index=False, sep='\t')
for level in set(calc_stats_df['level']):
	level_df = calc_stats_df.loc[calc_stats_df['level'] == level]
	sns.set_context("paper")
	ax = sns.catplot(x="statistic", y="score", hue='algorithm', kind='box',
						data=level_df, aspect=2
						)

	plt.plot([-1, 6], [0.25, 0.25], linestyle='--', alpha=0.3, color='k')
	plt.plot([-1, 6], [0.50, 0.50], linestyle='--', alpha=0.3, color='k')
	plt.plot([-1, 6], [0.75, 0.75], linestyle='--', alpha=0.3, color='k')

	plt.ylim(0, 1)
	plt.title('SAG-plus CAMI-1-High error analysis')
	ax._legend.set_title('Filter Type')
	plt.savefig(err_path + '/' + level + '_error_boxplox.svg',
				bbox_inches='tight'
				)
	plt.clf()

# build multi-level precision boxplot
level_list = ['genus', 'species', 'strain', 'perfect']
stat_list = ['precision', 'sensitivity', 'F1_score']
comb_stat_df = calc_stats_df.loc[((calc_stats_df['algorithm'] == 'combined') & 
									(calc_stats_df['level'].isin(level_list)) &
									(calc_stats_df['statistic'].isin(stat_list))
									)]
sns.set_context("paper")
ax = sns.catplot( x="level", y="score", hue='statistic', kind='box',
					data=comb_stat_df, aspect=2
					)

plt.plot([-1, 4], [0.25, 0.25], linestyle='--', alpha=0.3, color='k')
plt.plot([-1, 4], [0.50, 0.50], linestyle='--', alpha=0.3, color='k')
plt.plot([-1, 4], [0.75, 0.75], linestyle='--', alpha=0.3, color='k')

plt.ylim(0, 1)
plt.title('SAG-plus CAMI-1-High')

plt.savefig(err_path + '/multi-level_precision_boxplox.svg',
				bbox_inches='tight'
				)
plt.clf()
