import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import scatter
import math
from matplotlib.legend import Legend
import sys


working_dir = '/home/rmclaughlin/Ryan/Lulu/BinMulti/BM_190826/'
checkm_file = working_dir + 'MetaBAT2_BM_out_min1500_checkM_stdout_ALL.tsv'
checkm_df = pd.read_csv(checkm_file, sep='\t', header=0)
rpkm_file = working_dir + 'BM_out_binned.rpkm-1.csv'
rpkm_df = pd.read_csv(rpkm_file, sep=',', header=0)


print(checkm_df.head())
print(checkm_df.tail())

rpkm_df['bin'] = [x[1] + '.' + x[0].split('_')[2] if '_' in x[0]
					else x[1] + '.' + x[0] for x in
					zip(rpkm_df['Sequence_name'], rpkm_df['Sample'])
					]
rpkm_df['Sample'] = [x.split('.')[0] for x in rpkm_df['Sample']]
# Reproduce each bin for the great 8 (since they were binned as a co-asm)
wwtp_id_list = [x for x in set(rpkm_df['Sample']) if 'wastewater' in x]
wwtp_row_list = []
for i, row in checkm_df.loc[checkm_df['Sample'] == 'Great8'].iterrows():
	for w_id in wwtp_id_list:
		row['Sample'] = w_id
		wwtp_row_list.append(row)
# TODO: add wwtp to check_df and remove old great8 entries
checkm_df['CheckM Taxonomy'] = [x.rsplit(' ', 1)[0] for x in checkm_df['Marker lineage']]
checkm_df['bin'] = [x[1] + '.' + x[0].split('.')[1] for x in
					zip(checkm_df['Bin Id'], checkm_df['Sample'])
					]

group_rpkm_df = rpkm_df.groupby(['Sample', 'bin'])['RPKM'].sum().reset_index()
piv_rpkm_df = group_rpkm_df.pivot(index='bin', columns='Sample', values='RPKM').reset_index()
print(piv_rpkm_df.head())


merge_df = pd.merge(checkm_df, piv_rpkm_df, on='bin')
merge_df.to_csv(working_dir + 'lulu_bins_rpkms.tsv', sep='\t', header=True, index=False)

# build scatter of ALL
col_list = ['Bin Id', 'Marker lineage', '# genomes', '# markers', '# marker sets', '0',
				'1', '2', '3', '4', '5+', 'Completeness', 'Contamination',
				'Strain heterogeneity', 'Sample', 'CheckM Taxonomy'
				]
ex_col_list = [x for x in merge_df.columns if x not in col_list]
mean_rpkm_df = merge_df[ex_col_list].set_index('bin')
mean_rpkm_list = []
for i, row in mean_rpkm_df.iterrows():
	ave_r = row.sum()/row.count()
	mean_rpkm_list.append(ave_r)
merge_df['RPKM'] = mean_rpkm_list

RPKM_min = merge_df['RPKM'].min()
RPKM_max = merge_df['RPKM'].max()

magni = 10**(len(str(math.ceil(RPKM_max)))-1)
ru_max = int(math.ceil(RPKM_max/magni))* magni

tax_list = sorted(list(set(merge_df['CheckM Taxonomy'])), reverse=True)
# custom tax order
tax_list = ['root', 'k__Bacteria', 'k__Archaea', 'p__Proteobacteria', 'p__Firmicutes',
			'p__Euryarchaeota', 'p__Bacteroidetes', 'p__Actinobacteria', 'c__Spirochaetia',
			'c__Gammaproteobacteria', 'c__Deltaproteobacteria', 'c__Clostridia',
			'c__Betaproteobacteria', 'c__Alphaproteobacteria', 'o__Thermoanaerobacterales',
			'o__Sphingomonadales', 'o__Selenomonadales', 'o__Rickettsiales',
			'o__Rhizobiales', 'o__Lactobacillales', 'o__Flavobacteriales',
			'o__Cytophagales', 'o__Clostridiales', 'o__Burkholderiales',
			'o__Bacteroidales', 'o__Actinomycetales',  'f__Xanthomonadaceae',
			'f__Spirochaetaceae', 'f__Rhodocyclaceae', 'f__Rhodobacteraceae',
			'f__Moraxellaceae', 'f__Lachnospiraceae', 'f__Flavobacteriaceae',
			'f__Comamonadaceae'
			]
palette_list = sns.color_palette("Paired", n_colors=len(tax_list))
tax_pal_list = list(zip(tax_list, palette_list))
hue_dict = {x[0]:x[1] for x in tax_pal_list} 
sns.set_style("white")
sns.set_style("ticks")
sns.set_context("paper")
ax = sns.scatterplot(x='Completeness', y='Contamination', hue='CheckM Taxonomy',
						size='RPKM', edgecolor='gray', sizes=(0, ru_max/5), # sizes=(min_floor/25, max_ceil/25),
						data=merge_df, palette=hue_dict, alpha=0.75)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,  frameon=False)

legend_markers = [	
					scatter([-10], [0], marker='o',label=str(int(ru_max/100)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/20)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/10)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/4)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/2)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max)), color='k')
					]
lgnd = plt.legend(title='RPKM', handles=legend_markers, bbox_to_anchor=(1.6, 1), loc=2,
						borderaxespad=0., scatterpoints=1, fontsize=10, labelspacing=5,
						borderpad=3, handletextpad=3
						)
lgnd.legendHandles[0]._sizes = [int(int((ru_max/5)*0.01))]
lgnd.legendHandles[1]._sizes = [int(int((ru_max/5)*0.05))]
lgnd.legendHandles[2]._sizes = [int(int((ru_max/5)*0.1))]
lgnd.legendHandles[3]._sizes = [int(int((ru_max/5)*0.25))]
lgnd.legendHandles[4]._sizes = [int(int((ru_max/5)*0.50))]
lgnd.legendHandles[5]._sizes = [int(int(ru_max/5))]
plt.gca().add_artist(lgnd)

taxleg_markers = []
for t, p in tax_pal_list:
	s = scatter([-10], [0], marker='o',label=t, color=p, edgecolor='gray')
	taxleg_markers.append(s)

leg = plt.legend(title='CheckM_taxonomy', handles=taxleg_markers, bbox_to_anchor=(1.6, 1), loc=1,
					borderaxespad=0., scatterpoints=1, fontsize=10, labelspacing=1.25,
					borderpad=1
					)

for i, tm in enumerate(taxleg_markers):
	leg.legendHandles[i]._sizes = [75]

plt.gca().add_artist(leg)
plt.xlim(-5, 115)
plt.ylim(-30, 1000)
plt.savefig(working_dir + 'Lulu_Comp_Cont_ALL.svg', bbox_inches = 'tight')
plt.clf()

# Build scatter for MQ only
MQ_df = merge_df[(merge_df['Completeness'] >= 50) &
					(merge_df['Contamination'] <= 10)]

RPKM_min = MQ_df['RPKM'].min()
RPKM_max = MQ_df['RPKM'].max()
magni = 10**(len(str(math.ceil(RPKM_max)))-1)
ru_max = int(math.ceil(RPKM_max/magni))* magni

tax_list = sorted(list(set(MQ_df['CheckM Taxonomy'])), reverse=True)
#palette_list = sns.color_palette("Paired", n_colors=len(tax_list))
tax2_pal_list = [(x[0], x[1]) for x in tax_pal_list if x[0] in tax_list]
#hue_dict = {x[0]:x[1] for x in tax_pal_list} 
sns.set_style("white")
sns.set_style("ticks")
sns.set_context("notebook")
ax = sns.scatterplot(x='Completeness', y='Contamination', hue='CheckM Taxonomy',
						size='RPKM', edgecolor='gray', sizes=(0, ru_max/5), # sizes=(min_floor/25, max_ceil/25),
						data=MQ_df, palette=hue_dict, alpha=0.75)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,  frameon=False)

legend_markers = [	
					scatter([-10], [0], marker='o',label=str(int(ru_max/100)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/20)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/10)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/4)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/2)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max)), color='k')
					]
lgnd = plt.legend(title='RPKM', handles=legend_markers, bbox_to_anchor=(1.6, 1), loc=2,
						borderaxespad=0., scatterpoints=1, fontsize=10, labelspacing=5,
						borderpad=3, handletextpad=3
						)
lgnd.legendHandles[0]._sizes = [int(int((ru_max/5)*0.01))]
lgnd.legendHandles[1]._sizes = [int(int((ru_max/5)*0.05))]
lgnd.legendHandles[2]._sizes = [int(int((ru_max/5)*0.1))]
lgnd.legendHandles[3]._sizes = [int(int((ru_max/5)*0.25))]
lgnd.legendHandles[4]._sizes = [int(int((ru_max/5)*0.50))]
lgnd.legendHandles[5]._sizes = [int(int(ru_max/5))]
plt.gca().add_artist(lgnd)

taxleg_markers = []
for t, p in tax2_pal_list:
	s = scatter([-10], [0], marker='o',label=t, color=p, edgecolor='gray')
	taxleg_markers.append(s)

leg = plt.legend(title='CheckM_taxonomy', handles=taxleg_markers, bbox_to_anchor=(1.6, 1), loc=1,
					borderaxespad=0., scatterpoints=1, fontsize=10, labelspacing=1.25,
					borderpad=1
					)
for i, tm in enumerate(taxleg_markers):
	leg.legendHandles[i]._sizes = [75]

plt.gca().add_artist(leg)
plt.xlim(49, 101)
plt.ylim(-0.5, 11)
plt.savefig(working_dir + 'Lulu_Comp_Cont_MQHQ.svg', bbox_inches = 'tight')
plt.clf()


# Build scatter for Low Contamination only
LC_df = merge_df[(merge_df['Contamination'] <= 10)]

RPKM_min = LC_df['RPKM'].min()
RPKM_max = LC_df['RPKM'].max()
magni = 10**(len(str(math.ceil(RPKM_max)))-1)
ru_max = int(math.ceil(RPKM_max/magni))* magni

tax_list = sorted(list(set(LC_df['CheckM Taxonomy'])), reverse=True)
#palette_list = sns.color_palette("Paired", n_colors=len(tax_list))
tax2_pal_list = [(x[0], x[1]) for x in tax_pal_list if x[0] in tax_list]
#hue_dict = {x[0]:x[1] for x in tax_pal_list} 
sns.set_style("white")
sns.set_style("ticks")
sns.set_context("notebook")
ax = sns.scatterplot(x='Completeness', y='Contamination', hue='CheckM Taxonomy',
						size='RPKM', edgecolor='gray', sizes=(0, ru_max/5), # sizes=(min_floor/25, max_ceil/25),
						data=LC_df, palette=hue_dict, alpha=0.75)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,  frameon=False)

legend_markers = [	
					scatter([-10], [0], marker='o',label=str(int(ru_max/100)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/20)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/10)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/4)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/2)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max)), color='k')
					]
lgnd = plt.legend(title='RPKM', handles=legend_markers, bbox_to_anchor=(1.6, 1), loc=2,
						borderaxespad=0., scatterpoints=1, fontsize=10, labelspacing=5,
						borderpad=3, handletextpad=3
						)
lgnd.legendHandles[0]._sizes = [int(int((ru_max/5)*0.01))]
lgnd.legendHandles[1]._sizes = [int(int((ru_max/5)*0.05))]
lgnd.legendHandles[2]._sizes = [int(int((ru_max/5)*0.1))]
lgnd.legendHandles[3]._sizes = [int(int((ru_max/5)*0.25))]
lgnd.legendHandles[4]._sizes = [int(int((ru_max/5)*0.50))]
lgnd.legendHandles[5]._sizes = [int(int(ru_max/5))]
plt.gca().add_artist(lgnd)

taxleg_markers = []
for t, p in tax2_pal_list:
	s = scatter([-10], [0], marker='o',label=t, color=p, edgecolor='gray')
	taxleg_markers.append(s)

leg = plt.legend(title='CheckM_taxonomy', handles=taxleg_markers, bbox_to_anchor=(1.6, 1), loc=1,
					borderaxespad=0., scatterpoints=1, fontsize=10, labelspacing=1.25,
					borderpad=1
					)
for i, tm in enumerate(taxleg_markers):
	leg.legendHandles[i]._sizes = [75]

plt.gca().add_artist(leg)
plt.xlim(-5, 105)
plt.ylim(-0.5, 11)
plt.savefig(working_dir + 'Lulu_Comp_Cont_LC.svg', bbox_inches = 'tight')
plt.clf()


