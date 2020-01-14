import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import umap
import sys
import pandas as pd
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
import scipy.spatial as spatial

matrix_path = '/home/rmclaughlin/Ryan/Lulu/sourmash/MAG_distances.csv'
matrix_df = pd.read_csv(matrix_path, sep=',', header=0)
matrix_df.columns = [x.split('/')[8].split('.', 1)[0] + '.' + x.split('/')[8].rsplit('.', 2)[1] 
						for x in matrix_df.columns
						]
matrix_df.set_index(matrix_df.columns, inplace=True)
trans_df = matrix_df.T
AD_cols = [x for x in trans_df.columns if 'A' in x]
AD_rows = [x for x in trans_df.index if 'A' in x]
AD_matrix_df = trans_df[list(AD_cols)].loc[trans_df.index.isin(list(AD_rows))]

checkM_path = '/home/rmclaughlin/Ryan/Lulu/BinMulti/BM_190826/MetaBAT2_BM_out_min1500_checkM_stdout_ALL.tsv'
checkM_df = pd.read_csv(checkM_path, sep='\t', header=0)
alt_col_dict = {'High': 'orange', 'Medium': 'blue', 'Partial': 'gray', 'Low': 'white'}
alt_col_list = []
for i, row in checkM_df.iterrows():
	cont = row['Contamination']
	comp = row['Completeness']
	if (cont <= 5.0) & (comp >= 90.0):
		alt_col_list.append('High')
	elif ((cont <= 5.0) & (comp >= 50.0)):
		alt_col_list.append('Medium')
	elif (cont <= 5.0) & (comp < 50.0):
		alt_col_list.append('Partial')
	else:
		alt_col_list.append('Low')
checkM_df['MAG Quality'] = alt_col_list
checkM_df['Sample_Bin_ID'] = [x[0] + '.' + x[1].split('.')[1] for x in
								zip(checkM_df['Sample'], checkM_df['Bin Id'])
								]
AD_checkM_df = checkM_df.loc[checkM_df['Sample_Bin_ID'].isin(list(AD_matrix_df.index.values))]
HQ_checkM_df = AD_checkM_df.loc[AD_checkM_df['MAG Quality'] == 'High']


mag_map_dict = {x[0]:x[1] for x in zip(HQ_checkM_df['Sample_Bin_ID'],
					HQ_checkM_df['MAG Quality']) if x[1] not in
					['Low', 'Partial', 'Medium']
					}

MQ_cols = [x for x in AD_matrix_df.columns if x in mag_map_dict.keys()]
MQ_rows = [x for x in AD_matrix_df.index if x in mag_map_dict.keys()]
MQ_matrix_df = AD_matrix_df[list(MQ_cols)].loc[AD_matrix_df.index.isin(list(MQ_rows))]

grp_qual_list = [mag_map_dict[x] for x in AD_matrix_df.index.values
				if x in mag_map_dict.keys()
				]

#bincnt_dict = {x[0]: ('Singleton' if x[1] == 1 else 'Cluster') for x in
#				zip(MQ_matrix_df.index, MQ_matrix_df.ge(0.25).sum(axis=1))
#				}
bincnt_dict = {}
sm_cluster_lists = []
ge_df = MQ_matrix_df.ge(0.40)
for i in ge_df.index.values:
	sub_df = ge_df.loc[ge_df.index == i]
	sm_clus_list = [x[0] for x in zip(sub_df[sub_df].columns, sub_df.iloc[0])
							if x[1] == True
							]
	sm_cluster_lists.append(sm_clus_list)
	v_sum = len(sm_clus_list)
	key = i
	if v_sum == 1:
		val = '1'
	elif v_sum == 2:
		val = '2'
	elif v_sum == 3:
		val = '3'
	elif v_sum == 4:
		val = '4'
	elif v_sum > 4:
		val = '>4'
	bincnt_dict[key] = val
#sm_cluster_sets = set([tuple(x) for x in sm_cluster_lists])
#sm_clusters_delist = [sorted(list(x)) for x in sm_cluster_sets]
#sm_clusters_delist.sort(key=len, reverse=True)
delist_sets = {frozenset(e) for e in sm_cluster_lists}
new_delist = []
for e in delist_sets:
	if any(e < s for s in delist_sets):
		continue
	else:
		new_delist.append(list(e))

sm_cluster_dict = {}
c_index = 0
for clus in new_delist:
	for sam in clus:
		sm_cluster_dict[sam] = c_index
	c_index += 1

HQ_checkM_df['sm_clusters'] = [sm_cluster_dict[x] for x in HQ_checkM_df['Sample_Bin_ID']]

grp_sim_list = [bincnt_dict[x] for x in MQ_matrix_df.index.values
				if x in bincnt_dict.keys()
				]

MQ_matrix_df.to_csv('MAG_SM_Distance.tsv', sep='\t', index=True, header=True)
# TODO: mess with UMAP params to get best plot
sv_pth = './'
n_neighbors = 10
min_dist = 0.0
n_components = 2
metric = 'manhattan'
random_state = 42
spread = 5
title = ''
k_num = 50
fit = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, spread=spread,
				n_components=n_components, metric=metric, random_state=random_state
				)
features = MQ_matrix_df.values
targets = MQ_matrix_df.index.values

embedding = fit.fit_transform(features)
emb_df = pd.DataFrame(embedding, columns=['x', 'y'], index=MQ_matrix_df.index)
emb_df['BinID'] = [x for x in emb_df.index.values]
'''
emb_df = pd.DataFrame(embedding, columns=['x', 'y'], index=MQ_matrix_df.index).reset_index()
emb_df.columns = ['BinID', 'x', 'y']

# Group points by nearest neighbors
points = [tuple([x[1]['x'], x[1]['y']]) for x in emb_df.iterrows()]
point_tree = spatial.cKDTree(points)
cluster_dict = {}
clusters_list = list(point_tree.query_ball_point(points, 2).copy())
merged_cluster_lists = []
for ind, c_list in enumerate(clusters_list):
	merge_list = c_list.copy()
	p = False
	skip_list = []
	for ind2, c_list2 in enumerate(clusters_list):
		if ((ind != ind2) & (skip_list.count(ind) == 0)):
			for v in c_list:
				if ((c_list2.count(v) != 0) & (p == False)):
					merge_list.extend(c_list2)
					#skip_list.append(ind2)
					#p = True
	merged_cluster_lists.append(tuple(sorted(set(merge_list))))
'''

# plot by similarity
sns.despine()
sns.set(font='arial')
sns.set_style('white', {'axes.edgecolor': '1.0'})

fig, ax = plt.subplots()
# plot UMAP results
g = sns.scatterplot(x='x', y='y', data=emb_df, hue=grp_sim_list,
						hue_order=['1', '2', '3', '4', '>4']
						)
plt.gca().set_aspect('equal', 'datalim')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False)

#merged_cluster_set = list(set(merged_cluster_lists))
#count_list = [merged_cluster_set.count(x) for x in merged_cluster_set]
#for i, cluster in enumerate(merged_cluster_set):
for i, cluster in enumerate(new_delist):
	cluster_df = emb_df.loc[emb_df.index.isin(list(cluster))]
	#for binid in cluster_df['BinID']:
	#	cluster_dict[binid] = i
	data_df = cluster_df[['x', 'y']]
	data = [tuple([x[1]['x'], x[1]['y']]) for x in data_df.iterrows()]
	x, y = zip(*data)
	l = len(x)
	c_x = sum(x)/l
	c_y = sum(y)/l
	center = c_x, c_y
	radius = cdist(data_df, [center]).max() + 2
	ax.add_artist(plt.Circle(center, radius, zorder=0, alpha=0.25, edgecolor='k',
							facecolor='gray'))

#ax.axes.get_xaxis().set_visible(False)
#ax.axes.get_yaxis().set_visible(False)

plot_save_path = sv_pth + 'MAG_SM_Distance_single.pdf'
plt.savefig(plot_save_path, bbox_inches="tight")
plt.clf()


# plot by quality
sns.despine()
sns.set(font='arial')
sns.set_style('white', {'axes.edgecolor': '1.0'})

fig, ax = plt.subplots()
# plot UMAP results
g = sns.scatterplot(x='x', y='y', data=emb_df,
						hue=grp_qual_list
						)
plt.gca().set_aspect('equal', 'datalim')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False)


#merged_cluster_set = list(set(merged_cluster_lists))
#count_list = [merged_cluster_set.count(x) for x in merged_cluster_set]
#for i, cluster in enumerate(merged_cluster_set):
for i, cluster in enumerate(new_delist):
	cluster_df = emb_df.loc[emb_df.index.isin(list(cluster))]
	#for binid in cluster_df['BinID']:
	#	cluster_dict[binid] = i
	data_df = cluster_df[['x', 'y']]
	data = [tuple([x[1]['x'], x[1]['y']]) for x in data_df.iterrows()]
	x, y = zip(*data)
	l = len(x)
	c_x = sum(x)/l
	c_y = sum(y)/l
	center = c_x, c_y
	radius = cdist(data_df, [center]).max() + 2
	ax.add_artist(plt.Circle(center, radius, zorder=0, alpha=0.25, edgecolor='k',
							facecolor='gray'))
#emb_df['umap_clusters'] = [cluster_dict[x] for x in emb_df['BinID']]
emb_df['mem_cnt'] = [bincnt_dict[x] for x in emb_df['BinID']]
emb_df['qual_type'] = grp_qual_list


#ax.axes.get_xaxis().set_visible(False)
#ax.axes.get_yaxis().set_visible(False)

plot_save_path = sv_pth + 'MAG_SM_Distance_quality.pdf'
plt.savefig(plot_save_path, bbox_inches="tight")
plt.clf()
merge_df = pd.merge(emb_df, HQ_checkM_df, how='left', left_on='BinID',
					right_on='Sample_Bin_ID')
merge_df.sort_values(by=['sm_clusters'], ascending=True, inplace=True)
merge_df.to_csv('MAG_UMAP.tsv', sep='\t', index=False, header=True)
HQ_checkM_df['Completeness_int'] = [int(x) for x in HQ_checkM_df['Completeness']]

sorted_df = HQ_checkM_df.sort_values(by=['Completeness_int', 'Contamination'],
										ascending=[False, True]
										)
dedup_df = sorted_df.drop_duplicates(subset=['sm_clusters'], keep='first',
										inplace=False
										)
dedup_df.sort_values(by=['sm_clusters'], ascending=True, inplace=True)
dedup_df.to_csv('MAG_UMAP_reps.tsv', sep='\t', index=False, header=True)



'''
metaG_map_dict = {'11A': 'AD1', '13A': 'AD1', '15A': 'AD1', '15SI': 'WSS', '17A': 'AD1', '2AD43II': 'AD2',
					'2AD48II': 'AD2', '2AD52I': 'AD2','2AD57I': 'AD2', '2AD61I': 'AD2', '2AI': 'AD2',
					'2SI': 'WSS', '6A': 'AD1', '7A': 'AD1', '8A': 'AD1', '8SI': 'WSS', '9A': 'AD1',
					'9SII': 'WSS', 'AD118III': 'AD1', 'AD121II': 'AD1', 'AD126III': 'AD1',
					'AD128I': 'AD1', 'AD132III': 'AD1', 'AD138III': 'AD1', 'AD143II': 'AD1',
					'AD148III': 'AD1', 'AD152III': 'AD1', 'CEPT26II': 'CEPT',
					'SCT18I': 'WSS', 'SCT21II': 'WSS', 'SCT26II': 'WSS', 'SCT28I': 'WSS',
					'SCT32I': 'WSS', 'SCT38III': 'WSS', 'SCT43I': 'WSS', 'SCT48II': 'WSS',
					'SCT52I': 'WSS', 'SCT57I': 'WSS', 'SCT61I': 'WSS'
					}
grp_list = [metaG_map_dict[x] for x in filter2_df.index.values]
'''
'''
mag_map_dict = {'SCT57I_FD': 'WSS', '6A_II': 'AD1', 'AD152III_FD': 'AD1',
				'SCT61I_FD': 'WSS', 'CEPT26II_FD': 'CEPT', 'SCT28I_FD': 'WSS',
				'AD143II_FD': 'AD1', '15SI_FD': 'WSS', 'AD138III_FD': 'AD1',
				'2SI_FD': 'WSS', 'AD126III_FD': 'AD1', 'SCT38III_FD': 'WSS',
				'9SII_FD': 'WSS', 'AD121II_FD': 'AD1', 'AD118III_FD': 'AD1',
				'SCT52I_FD': 'WSS', '2AI_FD': 'AD2', 'SCT48II_FD': 'WSS', '9A_II': 'AD1',
				'SCT21II_FD': 'WSS', '2AD52I_FD': 'AD2', '11A_II': 'AD1', '2AD48II_FD': 'AD2',
				'2AD57I_FD': 'AD2', 'SCT26II_FD': 'WSS', 'SCT43I_FD': 'WSS', '2AD43II_FD': 'AD2',
				'AD148III_FD': 'AD1', '13A_III': 'AD1', '8A_II': 'AD1', 'AD128I_FD': 'AD1',
				'AD132III_FD': 'AD1', '2AD61I_FD': 'AD2', 'SCT18I_FD': 'WSS', '8SI_FD': 'WSS',
				'17A_III': 'AD1', '7A_III': 'AD1', '15A_II': 'AD1', 'SCT32I_FD': 'WSS'
				}
grp_list = [mag_map_dict[x.split('.')[0]] for x in matrix_df.index.values]
'''