import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import scatter
import math
from matplotlib.legend import Legend


checkm_file = 'MetaBAT2_lulu_min1500_checkM_stdout.tsv'
checkm_df = pd.read_csv(checkm_file, sep='\t', header=0)

rpkm_file = 'lulu_binned.rpkm.csv'
rpkm_df = pd.read_csv(rpkm_file, sep=',', header=0)


checkm_df['bin'] = [x.split('.')[1] for x in checkm_df['Bin Id']]
checkm_df['CheckM Taxonomy'] = [x.rsplit(' ', 1)[0] for x in checkm_df['Marker lineage']]
rpkm_df['bin'] = [x.split('_')[1] if '_' in x else x  for x in rpkm_df['Sequence_name']]
rpkm_df['Sample'] = [x.rsplit('_', 1)[0] for x in rpkm_df['Sample_name']]
group_rpkm_df = rpkm_df.groupby(['Sample', 'bin'])['RPKM'].sum().reset_index()
piv_rpkm_df = group_rpkm_df.pivot(index='bin', columns='Sample', values='RPKM').reset_index()

merge_df = pd.merge(checkm_df, piv_rpkm_df, on='bin')
merge_df.to_csv('lulu_bins_rpkms.tsv', sep='\t', header=True, index=False)

# build scatter of ALL
MQ_df = merge_df[(merge_df['Completeness'] >= 0) &
					(merge_df['Contamination'] <= 1000000)]
print(MQ_df.head())
MQ_df['RPKM'] = (MQ_df['wastewater_treat_10'] + MQ_df['wastewater_treat_11'] + MQ_df['wastewater_treat_12'] + \
				MQ_df['wastewater_treat_13'] + MQ_df['wastewater_treat_14'] + MQ_df['wastewater_treat_15'] + \
				MQ_df['wastewater_treat_16'] + MQ_df['wastewater_treat_9'])/8

RPKM_min = MQ_df['RPKM'].min()
RPKM_max = MQ_df['RPKM'].max()

magni = 10**(len(str(math.ceil(RPKM_max)))-1)
ru_max = int(math.ceil(RPKM_max/magni))* magni

tax_list = sorted(list(set(MQ_df['CheckM Taxonomy'])), reverse=True)
palette_list = sns.color_palette("Paired", n_colors=len(tax_list))
tax_pal_list = list(zip(tax_list, palette_list))
hue_dict = {x[0]:x[1] for x in tax_pal_list} 
sns.set_style("white")
sns.set_style("ticks")
sns.set_context("notebook")
ax = sns.scatterplot(x='Completeness', y='Contamination', hue='CheckM Taxonomy',
						size='RPKM', edgecolor='gray', sizes=(0, ru_max/5), # sizes=(min_floor/25, max_ceil/25),
						data=MQ_df, palette=hue_dict)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,  frameon=False)

legend_markers = [	
					scatter([-10], [0], marker='o',label=str(int(ru_max/100)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/20)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/10)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/4)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/2)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max)), color='k')
					]
lgnd = plt.legend(title='RPKM', handles=legend_markers, bbox_to_anchor=(1.5, 1), loc=2,
						borderaxespad=0., scatterpoints=1, fontsize=10, labelspacing=1.25,
						borderpad=1
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

leg = plt.legend(title='CheckM_taxonomy', handles=taxleg_markers, bbox_to_anchor=(1.5, 1), loc=1,
					borderaxespad=0., scatterpoints=1, fontsize=10, labelspacing=1.25,
					borderpad=1
					)

for i, tm in enumerate(taxleg_markers):
	leg.legendHandles[i]._sizes = [75]

plt.gca().add_artist(leg)
plt.xlim(-5, 105)
plt.ylim(-10, 600)
plt.savefig('Lulu_Comp_Cont_ALL.svg', bbox_inches = 'tight')
plt.clf()

# Build scatter for MQ only
MQ_df = merge_df[(merge_df['Completeness'] >= 50) &
					(merge_df['Contamination'] <= 10)]
print(MQ_df.head())
MQ_df['RPKM'] = (MQ_df['wastewater_treat_10'] + MQ_df['wastewater_treat_11'] + MQ_df['wastewater_treat_12'] + \
				MQ_df['wastewater_treat_13'] + MQ_df['wastewater_treat_14'] + MQ_df['wastewater_treat_15'] + \
				MQ_df['wastewater_treat_16'] + MQ_df['wastewater_treat_9'])/8

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
						data=MQ_df, palette=hue_dict)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,  frameon=False)

legend_markers = [	
					scatter([-10], [0], marker='o',label=str(int(ru_max/100)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/20)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/10)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/4)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max/2)), color='k'),
					scatter([-10], [0], marker='o',label=str(int(ru_max)), color='k')
					]
lgnd = plt.legend(title='RPKM', handles=legend_markers, bbox_to_anchor=(1.5, 1), loc=2,
						borderaxespad=0., scatterpoints=1, fontsize=10, labelspacing=1.25,
						borderpad=1
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

leg = plt.legend(title='CheckM_taxonomy', handles=taxleg_markers, bbox_to_anchor=(1.5, 1), loc=1,
					borderaxespad=0., scatterpoints=1, fontsize=10, labelspacing=1.25,
					borderpad=1
					)
for i, tm in enumerate(taxleg_markers):
	leg.legendHandles[i]._sizes = [75]

plt.gca().add_artist(leg)
plt.xlim(50, 100)
plt.ylim(-1, 11)
plt.savefig('Lulu_Comp_Cont_MQHQ.svg', bbox_inches = 'tight')
plt.clf()
