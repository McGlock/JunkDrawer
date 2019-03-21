import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import sys
from matplotlib.lines import Line2D
from matplotlib.pyplot import scatter
from sklearn import preprocessing
from operator import itemgetter



# Get collect all the LoopG data
looptax_path = '/home/rmclaughlin/LoopG/Ryan_work/sakinaw_work/loopG/clustering_analysis/quantification/'
looptax_prefix = 'all_depths.unique.good.sensitive_RC.unique.precluster.phylip.opti_mcc.0.03'
looptax_file_list = ['30M.cons.taxonomy', '36M.cons.taxonomy', '40M.cons.taxonomy',
						'45M.cons.taxonomy', '50M.cons.taxonomy', '60M.cons.taxonomy',
						'80M.cons.taxonomy', '120M.cons.taxonomy'
						]
looptax_df_list = []
for file in looptax_file_list:
	looptax_file = looptax_path + '.'.join([looptax_prefix, file])
	looptax_df = pd.read_csv(looptax_file, sep='\t', header=0)
	# remove singleton OTUs and calc proportions
	looptax_df = looptax_df.loc[looptax_df['Size'] > 1]
	#looptax_df['size_prop'] = (looptax_df['Size']/sum(looptax_df['Size']))*100
	#looptax_df = looptax_df.loc[looptax_df['size_prop'] > 0.1]
	# Add depth
	looptax_df['depth'] = file.split('.')[0]
	looptax_df_list.append(looptax_df)

loop_concat_df = pd.concat(looptax_df_list, ignore_index=True).reset_index()
#size = loop_concat_df[['Size']].values.astype(float)
#min_max_scaler = preprocessing.MinMaxScaler()
#size_scaled = min_max_scaler.fit_transform(size)
#loop_concat_df['Size'] = size_scaled
#print(loop_concat_df.groupby(['depth'])['Size'].sum())

# Get collect all the PacBio data
pacbiotax_path = '/home/rmclaughlin/LoopG/Ryan_work/sakinaw_work/pacbio/clustering_analysis/quantification/'
pacbiotax_prefix = 'all_depths.unique.good.sensitive_RC.unique.precluster.phylip.opti_mcc.0.03'
pacbiotax_file_list = ['30M.cons.taxonomy', '36M.cons.taxonomy', '40M.cons.taxonomy',
						'45M.cons.taxonomy', '50M.cons.taxonomy', '60M.cons.taxonomy',
						'80M.cons.taxonomy', '120M.cons.taxonomy'
						]
pacbiotax_df_list = []
for file in pacbiotax_file_list:
	pacbiotax_file = pacbiotax_path + '.'.join([pacbiotax_prefix, file])
	pacbiotax_df = pd.read_csv(pacbiotax_file, sep='\t', header=0)
	# remove singleton OTUs and calc proportions
	pacbiotax_df = pacbiotax_df.loc[pacbiotax_df['Size'] > 1]
	#pacbiotax_df['size_prop'] = (pacbiotax_df['Size']/sum(pacbiotax_df['Size']))*100
	#pacbiotax_df = pacbiotax_df.loc[pacbiotax_df['size_prop'] > 0.1]
	# Add depth
	pacbiotax_df['depth'] = file.split('.')[0]
	pacbiotax_df_list.append(pacbiotax_df)

pacbio_concat_df = pd.concat(pacbiotax_df_list, ignore_index=True).reset_index()
#print(pacbio_concat_df.groupby(['depth'])['Size'].sum())

# Add sequencing method to DFs
loop_concat_df['seqtype'] = 'LoopG'
pacbio_concat_df['seqtype'] = 'PacBio'
loop_concat_df.drop('index', axis=1, inplace=True)
pacbio_concat_df.drop('index', axis=1, inplace=True)

# Concat the two DFs
concat_df = pd.concat([loop_concat_df, pacbio_concat_df], ignore_index=True).reset_index()
concat_df.drop('index', axis=1, inplace=True)
taxonomy_list = [[y.rsplit('(', 1)[0] for y in x.split(';')]
					for x in list(set(concat_df['Taxonomy']))
					]

# pull in silva lineage into dataframe
silva_tax_df = pd.read_csv('/home/rmclaughlin/Hallam_Databases/formatted/mothur/' \
							'mothur/tax_slv_ssu_132.txt', sep='\t',
							names=['path', 'taxid', 'rank', 'remark', 'release']
							)
silva_tax_df['taxname'] = [x.rstrip(';').replace(' ', '_') if x.count(';') == 1
							else x.rstrip(';').rsplit(';', 1)[1].replace(' ', '_')
							for x in silva_tax_df['path']
							]
silva_tax_df.drop('path', axis=1, inplace=True)
taxpath_dict = {x:['domain', 'phylum', 'class', 'order', 'family', 'genus']
				for x in list(set(concat_df['Taxonomy']))
				}
for taxpath in taxpath_dict.keys():
	include_clade_list = taxpath_dict[taxpath]
	for i, clade in enumerate(include_clade_list):
		silva_clade_list = list(silva_tax_df.loc[
									silva_tax_df['rank'] == clade]['taxname']
									)
		found_clade_list = list(set([x for x in silva_clade_list if x in taxpath]))
		if (len(found_clade_list) > 1 and 'uncultured' in found_clade_list):
			found_clade_list.remove('uncultured')
			taxpath_dict[taxpath][i] = found_clade_list[0]
		elif len(found_clade_list) == 1:
			taxpath_dict[taxpath][i] = found_clade_list[0]
		elif len(found_clade_list) < 1:
			taxpath_dict[taxpath][i] = 'not_classified'

taxpath_df = pd.DataFrame.from_dict(taxpath_dict, orient='index',
									columns=['Domain','Phylum','Class',
									'Order', 'Family', 'Genus']
									)
taxpath_df.index.names = ['Taxonomy']
taxpath_df.reset_index(inplace=True)
merge_df = pd.merge(concat_df, taxpath_df, on='Taxonomy')
filter_df = merge_df[['OTU', 'Size', 'depth', 'seqtype', 'Domain', 'Phylum',
       'Class', 'Order', 'Family', 'Genus']]

'''
# Split the lineage info into major clades
splitter = lambda x: pd.Series([i.split('(')[0] for i in x.split(';')])
split_df = concat_df['Taxonomy'].apply(splitter)
split_df.rename(columns={0:'Domain', 1:'Phylum', 2:'Class', 3:'Order', 4:'Family',
						5:'Genus', 6:'Species'}, inplace=True)
# Drop species since its empty
#split_df.drop(columns=['Species'], inplace=True)

merge_df = pd.merge(concat_df[['OTU', 'Size', 'size_prop', 'seqtype', 'depth']],
						split_df, left_index=True, right_index=True
						)
# remove unclassified Bacteria
#merge_df = merge_df.loc[merge_df['Phylum'] != 'Bacteria_unclassified']
'''

# Order phylum by domain
dom_phy_list = zip(filter_df['Domain'], filter_df['Phylum'])
sort_dom_phy_list = sorted(set(dom_phy_list), key=itemgetter(0))
phylum_set = [x[1] for x in sort_dom_phy_list]

'''
# find all Candidate groups
for phylum in phylum_set:
	print(phylum)
	p_df = merge_df.loc[merge_df['Phylum'] == phylum]
	print(p_df.head())
	for col in p_df.columns:
		print(col)
		col_dat = p_df[col]
		print(list(col_dat.str.contains('andidat')))
'''

y_dict = dict((y[1], y[0]) for y in enumerate(phylum_set))
depth_order_list = ['30M', '36M', '40M', '45M', '50M', '60M', '80M', '120M']
x_dict = dict((x[1], x[0]) for x in enumerate(depth_order_list))

seqtype_list = list(set(filter_df['seqtype']))
domain_list = list(set(filter_df['Domain']))
for seqtype in seqtype_list:
	seqtype_df = filter_df.loc[filter_df['seqtype'] == seqtype]
	for domain in domain_list:
		domain_df = seqtype_df.loc[seqtype_df['seqtype'] == domain]
		domain_df['size_prop'] = (domain_df['Size']/sum(domain_df['Size']))*100

		phylum_df = domain_df.groupby(['Phylum', 'depth'])['size_prop'
					].sum().reset_index()

		phylum_df['x_index'] = [x_dict[x] for x in phylum_df['depth']]
		phylum_df['y_index'] = [y_dict[y] for y in phylum_df['Phylum']]

		legend_markers = [scatter([0], [0], marker='o',label='<0.1', color='k'),
						scatter([0], [0], marker='o',label='1', color='k'),
						scatter([0], [0], marker='o',label='5', color='k'),
						scatter([0], [0], marker='o',label='12.5', color='k'),
						scatter([0], [0], marker='o',label='25', color='k'),
						scatter([0], [0], marker='o',label='50', color='k')
						]
		plt.figure(figsize=(6,14))
		plt.grid(zorder=0)
		scaler = 8
		plt.scatter(x = phylum_df['x_index'],
		            y = phylum_df['y_index'],
		            s = phylum_df['size_prop']*scaler,
		            color='k',
		            alpha = 1, zorder=3)

		x_ticks = [x_dict[x] for x in x_dict.keys()]
		x_labels = [l for l in list(x_dict.keys())]
		plt.xticks(x_ticks, x_labels)

		y_ticks = [y_dict[y] for y in y_dict.keys()]
		y_labels = [l for l in list(y_dict.keys())]
		plt.yticks(y_ticks, y_labels)
		#plt.xticks(rotation=90)

		#plt.title(sample, fontsize=20)
		lgnd = plt.legend(title='Percent', handles=legend_markers, bbox_to_anchor=(1.05, 1),
							loc=2, borderaxespad=0., scatterpoints=1, fontsize=10,
							labelspacing=2, borderpad=4
							)
		lgnd.legendHandles[0]._sizes = [0.1*scaler]
		lgnd.legendHandles[1]._sizes = [1*scaler]
		lgnd.legendHandles[2]._sizes = [5*scaler]
		lgnd.legendHandles[3]._sizes = [12.5*scaler]
		lgnd.legendHandles[4]._sizes = [25*scaler]
		lgnd.legendHandles[5]._sizes = [50*scaler]
		#plt.xlabel('Sample ID')
		#plt.ylabel('Genus species')
		#plt.ylim(0, 100)
		plt.savefig('/home/rmclaughlin/LoopG/Ryan_work/code/' + seqtype +
						'_' + domain + 
						'_Depth_Phylum_Percent.svg', bbox_inches='tight'
						)
		plt.clf()

		phylum_df = domain_df.groupby(['Phylum', 'depth'])['Size'
					].sum().reset_index()

		phylum_df['x_index'] = [x_dict[x] for x in phylum_df['depth']]
		phylum_df['y_index'] = [y_dict[y] for y in phylum_df['Phylum']]

		legend_markers = [scatter([0], [0], marker='o',label='<10', color='k'),
						scatter([0], [0], marker='o',label='500', color='k'),
						scatter([0], [0], marker='o',label='1000', color='k'),
						scatter([0], [0], marker='o',label='2500', color='k'),
						scatter([0], [0], marker='o',label='5000', color='k'),
						scatter([0], [0], marker='o',label='10000', color='k')
						]
		plt.figure(figsize=(6,14))
		plt.grid(zorder=0)
		scaler = 0.05
		plt.scatter(x = phylum_df['x_index'],
		            y = phylum_df['y_index'],
		            s = phylum_df['Size']*scaler,
		            color='k',
		            alpha = 1, zorder=3)

		x_ticks = [x_dict[x] for x in x_dict.keys()]
		x_labels = [l for l in list(x_dict.keys())]
		plt.xticks(x_ticks, x_labels)

		y_ticks = [y_dict[y] for y in y_dict.keys()]
		y_labels = [l for l in list(y_dict.keys())]
		plt.yticks(y_ticks, y_labels)
		#plt.xticks(rotation=90)

		#plt.title(sample, fontsize=20)
		lgnd = plt.legend(title='# of OTUs', handles=legend_markers, bbox_to_anchor=(1.05, 1),
							loc=2, borderaxespad=0., scatterpoints=1, fontsize=10,
							labelspacing=2, borderpad=4
							)
		lgnd.legendHandles[0]._sizes = [10*scaler]
		lgnd.legendHandles[1]._sizes = [500*scaler]
		lgnd.legendHandles[2]._sizes = [1000*scaler]
		lgnd.legendHandles[3]._sizes = [2500*scaler]
		lgnd.legendHandles[4]._sizes = [5000*scaler]
		lgnd.legendHandles[5]._sizes = [10000*scaler]
		#plt.xlabel('Sample ID')
		#plt.ylabel('Genus species')
		#plt.ylim(0, 100)
		plt.savefig('/home/rmclaughlin/LoopG/Ryan_work/code/' +	seqtype +
					'_' + domain + 
					'_Depth_Phylum_raw_count.svg', bbox_inches='tight'
					)
		plt.clf()