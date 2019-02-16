import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import umap
import sys
import pandas as pd
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist


def plot_umap(df, grp_list, sv_pth='./', n_neighbors=15, min_dist=2,
				n_components=2, metric='euclidean', random_state=42,
				spread=2,
				title=''):
	fit = umap.UMAP(
		n_neighbors=n_neighbors,
		min_dist=min_dist,
		spread=spread,
		n_components=n_components,
		metric=metric,
		random_state=random_state
		)
	features = df.values
	targets = df.index.values

	embedding = fit.fit_transform(features)
	emb_df = pd.DataFrame(embedding, columns=['x', 'y'], index=df.index)
	# KMEANS for grouping on scatter
	kmeans = KMeans(4, random_state=0)
	labels = kmeans.fit(emb_df).predict(emb_df)
	plot_save_path = sv_pth + 'GroupBySampleType.svg'
	sns.despine()
	sns.set(font='arial')
	sns.set_style('white', {'axes.edgecolor': '1.0'})

	# plot UMAP results
	ax = sns.scatterplot(x='x', y='y', data=emb_df,
							hue=grp_list
							)
	plt.gca().set_aspect('equal', 'datalim')
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False)

	# plot the representation of the KMeans model
	centers = kmeans.cluster_centers_
	radii = [cdist(emb_df[labels == i], [center]).max()
				for i, center in enumerate(centers)]
	lp_dict = {0: (0,7), 1: (0,-8), 2: (0,5), 3: (0,5)}
	for i, v in enumerate(zip(centers, radii)):
		c, r = v[0], v[1]
		ax.add_patch(plt.Circle(c, r, zorder=0, alpha=0.25, edgecolor='k',
								facecolor='gray'))
		ax.annotate(i, xy=(c[0]+lp_dict[i][0], c[1]+lp_dict[i][1]),
						fontsize=12, ha="center"
						)
	ax.axes.get_xaxis().set_visible(False)
	ax.axes.get_yaxis().set_visible(False)
	plt.savefig(plot_save_path, bbox_inches="tight")
	plt.clf()

	emb_df['kmeans_groups'] = labels
	emb_df['Sample_type'] = grp_list
	emb_df.index.name = 'Sample_ID'
	emb_df.to_csv('GroupBySampleType.tsv', sep='\t', index=True, header=True)
	
	return emb_df


matrix = sys.argv[1]
metadata = sys.argv[2]

matrix_df = pd.read_csv(matrix, sep='\t', header=0)
otu_df = matrix_df.drop(['taxonomy'], axis=1).set_index('OTU ID')
#tax_df = matrix_df.drop(['OTU_ID'], axis=1)
metadata_df = pd.read_csv(metadata, sep='\t', header=0)
metadata_df.set_index('SampleID', inplace=True)
trans_otu_df = otu_df.T
joined_df = trans_otu_df.merge(metadata_df, how='outer', left_index=True, right_index=True)
ex_list = ['BarcodeSequence', 'LinkerPrimerSequence', 'UBC_ID', 'Sample_order',
			'Sample_type'
			]
group_list = joined_df['Sample_type']
features_df = joined_df[[x for x in joined_df.columns if x not in ex_list]]
plot_umap(features_df, group_list)
