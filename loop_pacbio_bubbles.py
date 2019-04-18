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
	#looptax_df = looptax_df.loc[looptax_df['Size'] > 1]
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
	#pacbiotax_df = pacbiotax_df.loc[pacbiotax_df['Size'] > 1]
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
taxpath_dict = {x:['domain', 'phylum', 'class', 'order',
					'family', 'genus'] for x in
					list(set(concat_df['Taxonomy']))
					}

for taxpath in taxpath_dict.keys():
	include_clade_list = taxpath_dict[taxpath]
	for i, clade in enumerate(include_clade_list):
		silva_clade_list = [x for x in list(silva_tax_df.loc[
									silva_tax_df['rank'] == clade]['taxname'])
									if x not in ['uncultured', 'unclassified']
									]
		taxpath_list = [y.rsplit('(', 1)[0] for y in taxpath.split(';')]
		found_clade_list = list(set([x for x in taxpath_list if x in silva_clade_list]))
		if len(found_clade_list) > 1:
			t_cnt_list = []
			for t in found_clade_list:
				silva_cnt = silva_clade_list.count(t)
				path_cnt = taxpath_list.count(t)
				if silva_cnt == path_cnt:
					new_list = [t]
			found_clade_list = new_list
			taxpath_dict[taxpath][i] = found_clade_list[0]
		elif len(found_clade_list) == 1:
			taxpath_dict[taxpath][i] = found_clade_list[0]
		elif len(found_clade_list) < 1:
			if i == 0:
				unclass_label = 'unclassified'
			elif i > 0:
				if (('_unclassified' in taxpath_dict[taxpath][i-1]) or
					('unclassified' in taxpath_dict[taxpath][i-1])
					):
					unclass_label = taxpath_dict[taxpath][i-1]
				else: 
					unclass_label = taxpath_dict[taxpath][i-1] + '_unclassified'
			taxpath_dict[taxpath][i] = unclass_label

taxpath_df = pd.DataFrame.from_dict(taxpath_dict, orient='index',
									columns=['Domain','Phylum','Class',
									'Order', 'Family', 'Genus']
									)
taxpath_df.index.names = ['Taxonomy']
taxpath_df.reset_index(inplace=True)
merge_df = pd.merge(concat_df, taxpath_df, on='Taxonomy')
filter_df = merge_df[['OTU', 'Size', 'depth', 'seqtype', 'Domain', 'Phylum',
       'Class', 'Order', 'Family', 'Genus'
       ]]

# find all Candidate groups
candidate_key_list = ['Gracilibacteria', 'Aerophobetes', 'BRC1', 'Fermentibacteria',
						'JL-ETNP-Z39', 'Kazan-3B-28', 'LD1-PA38', 'MVP-21', 'NPL-UPA2',
						'OC31', 'Parcubacteria', 'Microgenomates', 'Omnitrophica',
						'Aminicenantes', 'Atribacteria', 'RF3', 'SM2F11',
						'Absconditabacteria', 'TA06', 'Dependentiae', 'Saccharibacteria',
						'WCHB1-60', 'Latescibacteria', 'WS6', 'Cloacimonetes', 'BD1-5',
						'CD12', 'BRC1', 'Hyd24-12', 'JL-ETNP-Z39', 'Kazan-3B-28',
						'LD1-PA38', 'MVP-21', 'NPL-UPA2', 'OC31', 'OD1', 'OP11',
						'OP3', 'OP8', 'OP9/JS1', 'RF3', 'SM2F11', 'SR1', 'TA06',
						'TM6', 'TM7', 'WCHB1-60', 'WS3', 'WS6', 'WWE1'
						]
candidate_list = []
for i, row in filter_df.iterrows():
	inter_list = [x for x in row if str(x) in candidate_key_list]
	if ((any('andidat' in str(s) for s in row)) or(len(inter_list) != 0)):
		candidate_list.append(True)
	else:
		candidate_list.append(False)
filter_df['candidate'] = candidate_list
# join phylum and class
filter_df['Phy_Cla'] = [x[0] + ' ' + x[1] for x in
						zip(filter_df['Phylum'], filter_df['Class'])
						]
# add * to all candidate phyla
filter_df['Phyl_Cla_label'] = ['*'+x[0]+'*' if x[1] == True else x[0]
						for x in zip(filter_df['Phy_Cla'], filter_df['candidate'])
						]
filter_df['Phyl_label'] = ['*'+x[0]+'*' if x[1] == True else x[0]
						for x in zip(filter_df['Phylum'], filter_df['candidate'])
						]

# remove the Euks
filter_df = filter_df[filter_df['Domain'] != 'Eukaryota']

# Order phylum by domain
dom_phy_list = zip(filter_df['Domain'], filter_df['Phyl_label'])
sort_dom_phy_list = sorted(set(dom_phy_list), key=itemgetter(1), reverse=True)
sort_dom_phy_list.sort(key=itemgetter(0))
phylum_set = [x[1] for x in sort_dom_phy_list if '*' not in x[1]]
# add Candidates to end
#cand_phy_list = [x[1] for x in sort_dom_phy_list if '*' in x[1]]
#cand_phy_list.extend(phylum_set)
#phylum_set = cand_phy_list

# Order class by domain
dom_cla_list = zip(filter_df['Domain'], filter_df['Phyl_Cla_label'])
sort_dom_cla_list = sorted(set(dom_cla_list), key=itemgetter(1), reverse=True)
sort_dom_cla_list.sort(key=itemgetter(0))
class_set = [x[1] for x in sort_dom_cla_list if '*' not in x[1]]
# add Candidates to end
#cand_cla_list = [x[1] for x in sort_dom_cla_list if '*' in x[1]]
#cand_cla_list.extend(class_set)
#class_set = cand_cla_list

seqtype_list = list(set(filter_df['seqtype']))
domain_list = list(set(filter_df['Domain']))
for seqtype in seqtype_list:
	print(seqtype)
	seqtype_df = filter_df.loc[(filter_df['seqtype'] == seqtype)]
	#domain_df = seqtype_df.copy()
	#for domain in domain_list:
		#domain_df = seqtype_df.loc[seqtype_df['Domain'] == domain]
		#if domain_df.empty != True:
	depth_sum_dict = seqtype_df.groupby(['depth'])['Size'].sum().to_dict()
	seqtype_df['depth_sum'] = [depth_sum_dict[x] for x in seqtype_df['depth']]
	seqtype_df['size_prop'] = (seqtype_df['Size']/seqtype_df['depth_sum'])*100
	#seqtype_df = seqtype_df.loc[seqtype_df['Size'] > 1]

	# group candiate groups
	cand_phy_df = seqtype_df.loc[seqtype_df['candidate'] == True]
	grp_cand_df = cand_phy_df.groupby(['depth'])['size_prop', 'Size'].sum().reset_index()
	grp_cand_df['Phylum'] = 'Candidate Groups'
	grp_cand_df['Phyl_label'] = 'Candidate Groups'
	grp_cand_df['candidate'] = True
	phylum_set.insert(0, 'Candidate Groups')
	class_set.insert(0, 'Candidate Groups')

	# add grouped candidates to final df
	add_cand_df = pd.concat([seqtype_df.loc[seqtype_df['candidate'] == False], grp_cand_df])

	# build class and phylum level df	
	class_df = add_cand_df.groupby(['Phylum', 'Phyl_Cla_label', 'depth', 'candidate']
									)['size_prop', 'Size'].sum().reset_index()
	phylum_df = add_cand_df.groupby(['Phylum', 'Phyl_label', 'depth', 'candidate']
									)['size_prop', 'Size'].sum().reset_index()			

	phylum_set.insert(0, '')
	phylum_set.append('')
	class_set.insert(0, '')
	class_set.append('')

	# Phylum % plot
	phy_y_dict = dict((y[1], y[0]) for y in enumerate(phylum_set))
	depth_order_list = ['30M', '36M', '40M', '45M', '50M', '60M', '80M', '120M']
	x_dict = dict((x[1], x[0]) for x in enumerate(depth_order_list))

	phylum_df['x_index'] = [x_dict[x] for x in phylum_df['depth']]
	phylum_df['y_index'] = [phy_y_dict[y] for y in phylum_df['Phyl_label']]

	legend_markers = [scatter([0], [0], marker='o',label='<0.1', color='k'),
					scatter([0], [0], marker='o',label='1', color='k'),
					scatter([0], [0], marker='o',label='10', color='k'),
					scatter([0], [0], marker='o',label='25', color='k'),
					scatter([0], [0], marker='o',label='50', color='k'),
					scatter([0], [0], marker='o',label='100', color='k')
					]
	plt.figure(figsize=(6,18))
	plt.grid(zorder=0)
	scaler = 5
	plt.scatter(x = phylum_df['x_index'],
	            y = phylum_df['y_index'],
	            s = phylum_df['size_prop']*scaler,
	            color='k',
	            alpha = 1, zorder=3)

	x_ticks = [x_dict[x] for x in x_dict.keys()]
	x_labels = [l for l in list(x_dict.keys())]
	plt.xticks(x_ticks, x_labels)

	y_ticks = [phy_y_dict[y] for y in phy_y_dict.keys()]
	y_labels = [l for l in list(phy_y_dict.keys())]
	plt.tick_params(axis='both', labelsize=10)
	plt.yticks(y_ticks, y_labels)
	lgnd = plt.legend(title='Percent', handles=legend_markers, bbox_to_anchor=(1.05, 1),
						loc=2, borderaxespad=0., scatterpoints=1, fontsize=10,
						labelspacing=2, borderpad=4
						)
	lgnd.legendHandles[0]._sizes = [0.1*scaler]
	lgnd.legendHandles[1]._sizes = [1*scaler]
	lgnd.legendHandles[2]._sizes = [10*scaler]
	lgnd.legendHandles[3]._sizes = [25*scaler]
	lgnd.legendHandles[4]._sizes = [50*scaler]
	lgnd.legendHandles[5]._sizes = [100*scaler]
	plt.savefig('/home/rmclaughlin/LoopG/Ryan_work/code/' + seqtype +
					'_Depth_Phylum_Percent.svg', bbox_inches='tight'
					)
	plt.clf()

	# Class % plot
	cla_y_dict = dict((y[1], y[0]) for y in enumerate(class_set))
	depth_order_list = ['30M', '36M', '40M', '45M', '50M', '60M', '80M', '120M']
	x_dict = dict((x[1], x[0]) for x in enumerate(depth_order_list))

	class_df['x_index'] = [x_dict[x] for x in class_df['depth']]
	class_df['y_index'] = [cla_y_dict[y] for y in class_df['Phyl_Cla_label']]

	legend_markers = [scatter([0], [0], marker='o',label='<0.1', color='k'),
					scatter([0], [0], marker='o',label='1', color='k'),
					scatter([0], [0], marker='o',label='10', color='k'),
					scatter([0], [0], marker='o',label='25', color='k'),
					scatter([0], [0], marker='o',label='50', color='k'),
					scatter([0], [0], marker='o',label='100', color='k')
					]
	plt.figure(figsize=(6,24))
	plt.grid(zorder=0)
	scaler = 5
	plt.scatter(x = class_df['x_index'],
	            y = class_df['y_index'],
	            s = class_df['size_prop']*scaler,
	            color='k',
	            alpha = 1, zorder=3)

	x_ticks = [x_dict[x] for x in x_dict.keys()]
	x_labels = [l for l in list(x_dict.keys())]
	plt.xticks(x_ticks, x_labels)

	y_ticks = [cla_y_dict[y] for y in cla_y_dict.keys()]
	y_labels = [l for l in list(cla_y_dict.keys())]
	plt.tick_params(axis='both', labelsize=10)
	plt.yticks(y_ticks, y_labels)
	lgnd = plt.legend(title='Percent', handles=legend_markers, bbox_to_anchor=(1.05, 1),
						loc=2, borderaxespad=0., scatterpoints=1, fontsize=10,
						labelspacing=2, borderpad=4
						)
	lgnd.legendHandles[0]._sizes = [0.1*scaler]
	lgnd.legendHandles[1]._sizes = [1*scaler]
	lgnd.legendHandles[2]._sizes = [10*scaler]
	lgnd.legendHandles[3]._sizes = [25*scaler]
	lgnd.legendHandles[4]._sizes = [50*scaler]
	lgnd.legendHandles[5]._sizes = [100*scaler]
	plt.savefig('/home/rmclaughlin/LoopG/Ryan_work/code/' + seqtype +
					'_Depth_Class_Percent.svg', bbox_inches='tight'
					)
	plt.clf()

	# Phylum Raw Count
	legend_markers = [scatter([0], [0], marker='o',label='<10', color='k'),
					scatter([0], [0], marker='o',label='500', color='k'),
					scatter([0], [0], marker='o',label='1000', color='k'),
					scatter([0], [0], marker='o',label='2500', color='k'),
					scatter([0], [0], marker='o',label='5000', color='k'),
					scatter([0], [0], marker='o',label='10000', color='k')
					]
	plt.figure(figsize=(6,18))
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

	y_ticks = [phy_y_dict[y] for y in phy_y_dict.keys()]
	y_labels = [l for l in list(phy_y_dict.keys())]
	plt.tick_params(axis='both', labelsize=10)
	plt.yticks(y_ticks, y_labels)
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
	plt.savefig('/home/rmclaughlin/LoopG/Ryan_work/code/' +	seqtype +
				'_Depth_Phylum_raw_count.svg', bbox_inches='tight'
				)
	plt.clf()


	# Class Raw Count
	legend_markers = [scatter([0], [0], marker='o',label='<10', color='k'),
					scatter([0], [0], marker='o',label='500', color='k'),
					scatter([0], [0], marker='o',label='1000', color='k'),
					scatter([0], [0], marker='o',label='2500', color='k'),
					scatter([0], [0], marker='o',label='5000', color='k'),
					scatter([0], [0], marker='o',label='10000', color='k')
					]
	plt.figure(figsize=(6,24))
	plt.grid(zorder=0)
	scaler = 0.05
	plt.scatter(x = class_df['x_index'],
	            y = class_df['y_index'],
	            s = class_df['Size']*scaler,
	            color='k',
	            alpha = 1, zorder=3)

	x_ticks = [x_dict[x] for x in x_dict.keys()]
	x_labels = [l for l in list(x_dict.keys())]
	plt.xticks(x_ticks, x_labels)

	y_ticks = [cla_y_dict[y] for y in cla_y_dict.keys()]
	y_labels = [l for l in list(cla_y_dict.keys())]
	plt.tick_params(axis='both', labelsize=10)
	plt.yticks(y_ticks, y_labels)
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
	plt.savefig('/home/rmclaughlin/LoopG/Ryan_work/code/' +	seqtype +
				'_Depth_Class_raw_count.svg', bbox_inches='tight')
	plt.clf()
