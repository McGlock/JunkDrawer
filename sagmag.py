import sys
from os import listdir
from os.path import isfile, join
from itertools import islice
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import hdbscan
import umap
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.datasets import load_iris, load_digits
from sklearn.model_selection import train_test_split
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
import sklearn.cluster as cluster
from sklearn.mixture import GaussianMixture as GMM
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import LeaveOneOut
from sklearn.neighbors import KernelDensity
from sklearn.base import BaseEstimator, ClassifierMixin
import pickle



sns.set(style='white', context='notebook', rc={'figure.figsize':(14,10)})


def tetra_cnt(seq_list):
	tetra_cnt_dict = {'aaaa': [], 'aaat': [], 'aaag': [], 'aaac': [], 'aata': [],
					'aatt': [], 'aatg': [], 'aatc': [], 'aaga': [], 'aagt': [],
					'aagg': [], 'aagc': [], 'aaca': [], 'aact': [], 'aacg': [],
					'aacc': [], 'ataa': [], 'atat': [], 'atag': [], 'atac': [],
					'atta': [], 'attt': [], 'attg': [], 'attc': [], 'atga': [],
					'atgt': [], 'atgg': [], 'atgc': [], 'atca': [], 'atct': [], 
					'atcg': [], 'atcc': [], 'agaa': [], 'agat': [], 'agag': [],
					'agac': [], 'agta': [], 'agtt': [], 'agtg': [], 'agtc': [],
					'agga': [], 'aggt': [], 'aggg': [], 'aggc': [], 'agca': [],
					'agct': [], 'agcg': [], 'agcc': [], 'acaa': [], 'acat': [],
					'acag': [], 'acac': [], 'acta': [], 'actt': [], 'actg': [],
					'actc': [], 'acga': [], 'acgt': [], 'acgg': [], 'acgc': [], 
					'acca': [], 'acct': [], 'accg': [], 'accc': [], 'taaa': [],
					'taat': [], 'taag': [], 'taac': [], 'tata': [], 'tatt': [],
					'tatg': [], 'tatc': [], 'taga': [], 'tagt': [], 'tagg': [],
					'tagc': [], 'taca': [], 'tact': [], 'tacg': [], 'tacc': [],
					'ttaa': [], 'ttat': [], 'ttag': [], 'ttac': [], 'ttta': [],
					'tttt': [], 'tttg': [], 'tttc': [], 'ttga': [], 'ttgt': [], 
					'ttgg': [], 'ttgc': [], 'ttca': [], 'ttct': [], 'ttcg': [],
					'ttcc': [], 'tgaa': [], 'tgat': [], 'tgag': [], 'tgac': [],
					'tgta': [], 'tgtt': [], 'tgtg': [], 'tgtc': [], 'tgga': [],
					'tggt': [], 'tggg': [], 'tggc': [], 'tgca': [], 'tgct': [],
					'tgcg': [], 'tgcc': [], 'tcaa': [], 'tcat': [], 'tcag': [],
					'tcac': [], 'tcta': [], 'tctt': [], 'tctg': [], 'tctc': [], 
					'tcga': [], 'tcgt': [], 'tcgg': [], 'tcgc': [], 'tcca': [],
					'tcct': [], 'tccg': [], 'tccc': [], 'gaaa': [], 'gaat': [],
					'gaag': [], 'gaac': [], 'gata': [], 'gatt': [], 'gatg': [],
					'gatc': [], 'gaga': [], 'gagt': [], 'gagg': [], 'gagc': [],
					'gaca': [], 'gact': [], 'gacg': [], 'gacc': [], 'gtaa': [],
					'gtat': [], 'gtag': [], 'gtac': [], 'gtta': [], 'gttt': [], 
					'gttg': [], 'gttc': [], 'gtga': [], 'gtgt': [], 'gtgg': [],
					'gtgc': [], 'gtca': [], 'gtct': [], 'gtcg': [], 'gtcc': [],
					'ggaa': [], 'ggat': [], 'ggag': [], 'ggac': [], 'ggta': [],
					'ggtt': [], 'ggtg': [], 'ggtc': [], 'ggga': [], 'gggt': [],
					'gggg': [], 'gggc': [], 'ggca': [], 'ggct': [], 'ggcg': [],
					'ggcc': [], 'gcaa': [], 'gcat': [], 'gcag': [], 'gcac': [], 
					'gcta': [], 'gctt': [], 'gctg': [], 'gctc': [], 'gcga': [],
					'gcgt': [], 'gcgg': [], 'gcgc': [], 'gcca': [], 'gcct': [],
					'gccg': [], 'gccc': [], 'caaa': [], 'caat': [], 'caag': [],
					'caac': [], 'cata': [], 'catt': [], 'catg': [], 'catc': [],
					'caga': [], 'cagt': [], 'cagg': [], 'cagc': [], 'caca': [],
					'cact': [], 'cacg': [], 'cacc': [], 'ctaa': [], 'ctat': [], 
					'ctag': [], 'ctac': [], 'ctta': [], 'cttt': [], 'cttg': [],
					'cttc': [], 'ctga': [], 'ctgt': [], 'ctgg': [], 'ctgc': [],
					'ctca': [], 'ctct': [], 'ctcg': [], 'ctcc': [], 'cgaa': [],
					'cgat': [], 'cgag': [], 'cgac': [], 'cgta': [], 'cgtt': [],
					'cgtg': [], 'cgtc': [], 'cgga': [], 'cggt': [], 'cggg': [],
					'cggc': [], 'cgca': [], 'cgct': [], 'cgcg': [], 'cgcc': [], 
					'ccaa': [], 'ccat': [], 'ccag': [], 'ccac': [], 'ccta': [], 
					'cctt': [], 'cctg': [], 'cctc': [], 'ccga': [], 'ccgt': [],
					'ccgg': [], 'ccgc': [], 'ccca': [], 'ccct': [], 'cccg': [],
					'cccc': []
					}  # build empty dict or tetranucleotide counting

	# count up all kmers and also populate the tetra dict
	for seq in seq_list:
		tmp_dict = {k: 0 for k, v in tetra_cnt_dict.items()}
		total_kmer_cnt = 0
		clean_seq = seq.strip('\n').lower()
		kmer_list = [''.join(x) for x in get_kmer(clean_seq, 4)]
		for tetra in tmp_dict.keys():
			count_tetra = kmer_list.count(tetra)
			tmp_dict[tetra] = count_tetra
			total_kmer_cnt += count_tetra
		# map tetras to their reverse tetras (not compliment)
		dedup_dict = {}
		for tetra in tmp_dict.keys():
			if (tetra not in dedup_dict.keys()) & (tetra[::-1] not in dedup_dict.keys()):
				dedup_dict[tetra] = ''
			elif tetra[::-1] in dedup_dict.keys():
				dedup_dict[tetra[::-1]] = tetra
		# combine the tetras and their reverse (not compliment), convert to proportions
		tetra_prop_dict = {}
		for tetra in dedup_dict.keys():
			if dedup_dict[tetra] != '':
				#tetra_prop_dict[tetra] = tmp_dict[tetra] + tmp_dict[dedup_dict[tetra]]
				t_prop = (tmp_dict[tetra] 
							+ tmp_dict[dedup_dict[tetra]]) / total_kmer_cnt
				tetra_prop_dict[tetra] = t_prop
			else:
				#tetra_prop_dict[tetra] = tmp_dict[tetra]
				t_prop = tmp_dict[tetra] / total_kmer_cnt
				tetra_prop_dict[tetra] = t_prop
		# add to tetra_cnt_dict
		for k in tetra_cnt_dict.keys():
			if k in tetra_prop_dict.keys():
				tetra_cnt_dict[k].append(tetra_prop_dict[k])
			else:
				tetra_cnt_dict[k].append(0.0)
	# convert the final dict into a pd dataframe for ease
	tetra_cnt_df = pd.DataFrame.from_dict(tetra_cnt_dict)
	dedupped_df = tetra_cnt_df.loc[:, (tetra_cnt_df != 0.0).any(axis=0)]
	return dedupped_df


def get_kmer(seq, n):  # found on internet
	"Returns a sliding window (of width n) over data from the iterable"
	"   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...				   "
	it = iter(seq)
	result = tuple(islice(it, n))
	if len(result) == n:
		yield result
	for elem in it:
		result = result[1:] + (elem,)
		yield result


def get_frags(seq, l_max, o_lap):  # not sure about this function, probs shouldn't use it
	"Fragments seq into subseqs of length l_max and overlap of o_lap"
	"Leftover tail overlaps with tail-1"
	seq_frags = []
	if (l_max != 0) and (len(seq) > l_max):
		offset = l_max - o_lap
		for i in range(0, len(seq), offset):
			if i+l_max < len(seq):
				frag = seq[i:i+l_max]
			else:
				frag = seq[-l_max:]
			if 'n' not in frag:
				seq_frags.append(frag)
	elif 'n' not in seq:
		seq_frags.append(seq)
	else:
		print('You have Ns in your seqs :(')
	return seq_frags

def get_subseqs(seq_list, n, o_lap):
	all_sub_seqs = []
	all_sub_headers = []
	for seq_tup in seq_list:
		header, seq = seq_tup
		clean_seq = seq.strip('\n').lower()
		sub_list = get_frags(clean_seq, n, o_lap)
		sub_headers = [header + '_' + str(i) for i, x in enumerate(sub_list, start=0)]
		all_sub_seqs.extend(sub_list)
		all_sub_headers.extend(sub_headers)	
	return all_sub_headers, all_sub_seqs


def get_seqs(fasta_file):
	sag_contigs = []
	with open(fasta_file, 'r') as f:
		data = f.read()
		split_data = data.split('>')
		for record in split_data:
			split_rec = record.split('\n')
			header = split_rec[0]
			seq = ''.join(split_rec[1:])
			if seq != '':
				sag_contigs.append((header, seq))
	return sag_contigs


def calc_seg(subseqs):
	seg_list = []
	for seq in subseqs:
		seg_sum = 0
		for i, nuc in enumerate(seq, start=0):
			seg_sum += calc_nuc(nuc, i)
		seg_list.append(seg_sum)
	return seg_list


def calc_nuc(nuc, ind):
	nuc_dict = {'a': 0, 't': 1, 'c': 2, 'g': 3}
	nuc_hash = nuc_dict[nuc] * (4**ind)
	return nuc_hash


def plot_umap(df, sv_pth='./', n_neighbors=15, min_dist=0.1,
				n_components=2, metric='euclidean', random_state=42,
				title=''):
	fit = umap.UMAP(
		n_neighbors=n_neighbors,
		min_dist=min_dist,
		n_components=n_components,
		metric=metric,
		random_state=random_state
		)
	features = df.values
	targets = df.index.values
	targets_ints = [x[0] for x in enumerate(targets, start=0)]

	embedding = fit.fit_transform(features)
	ax = sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], hue=targets)
	plt.gca().set_aspect('equal', 'datalim')
	plt.title(title, fontsize=18)
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plot_file_name = '_'.join([str(n_neighbors), str(min_dist),
								str(n_components), metric]) + '.png'
	plot_save_path = join(sv_pth, plot_file_name)
	plt.savefig(plot_save_path, bbox_inches="tight")
	plt.clf()
	return embedding


def draw_ellipse(position, covariance, ax=None, **kwargs):
	"""Draw an ellipse with a given position and covariance"""
	ax = ax or plt.gca()
	
	# Convert covariance to principal axes
	if covariance.shape == (2, 2):
		U, s, Vt = np.linalg.svd(covariance)
		angle = np.degrees(np.arctan2(U[1, 0], U[0, 0]))
		width, height = 2 * np.sqrt(s)
	else:
		angle = 0
		width, height = 2 * np.sqrt(covariance)
	
	# Draw the Ellipse
	for nsig in range(1, 4):
		ax.add_patch(Ellipse(position, nsig * width, nsig * height,
							 angle, **kwargs))
		
def plot_gmm(sg_id, sv_pth, gmm, X, targets, label=True, ax=None):
	labels = gmm.fit(X).predict(X)
	if label:
		ax = sns.scatterplot(x=data[:, 0], y=data[:, 1], hue=targets)
	plt.gca().set_aspect('equal', 'datalim')
	
	w_factor = 0.2 / gmm.weights_.max()
	for pos, covar, w in zip(gmm.means_, gmm.covariances_, gmm.weights_):
		draw_ellipse(pos, covar, alpha=w * w_factor)
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plot_file_name = sag_id + '.' + 'GMM_predict.UMAP.png'
	plot_save_path = join(save_path, plot_file_name)
	plt.savefig(plot_save_path, bbox_inches="tight")
	plt.clf()



### Start Main ###
#sag_path = sys.argv[1]
sag_file = sys.argv[1]
mg_file = sys.argv[2]
max_contig_len = int(sys.argv[3])
overlap_len = int(sys.argv[4])
save_path = sys.argv[5]

if max_contig_len > 0:
	print('Max contig size is %s bp' % max_contig_len)
else:
	print('Not fragmenting contigs')

# SAG Tetras
with open(sag_file.split('.')[0] + '.full.headers', 'r') as f:
    sag_raw_contig_headers = [x.replace('>', '') for x in f.read().splitlines()]
sag_id = sag_file.split('/')[-1].split('.')[0]
# Break up contigs into overlapping subseqs
sag_contigs = get_seqs(sag_file)
sag_headers, sag_subs = get_subseqs(sag_contigs, max_contig_len, overlap_len)
sag_tetra_df = pd.DataFrame.from_dict(tetra_cnt(sag_subs))
sag_tetra_df['contig_id'] = sag_headers
sag_tetra_df.set_index('contig_id', inplace=True)
print('SAG %s tetranucleotide frequencies calculated' % sag_id)
sag_tetra_df.to_csv(join(save_path, sag_id + '.tsv'), sep='\t')
'''
sag_tetra_df = pd.read_csv(join(save_path, sag_id + '.tsv'), sep='\t', index_col=0,
						header=0)
print('Opened SAG tetranucleotide tsv file')
'''
# SAG L-seg hash
tmp, sag_Ls = get_subseqs(sag_contigs, 24, 23)
sag_hashes = calc_seg(sag_Ls)
sag_hashes.sort(reverse=True)
sag_hashes_set = set(sag_hashes)
with open(join(save_path, sag_id + '.pkl'), 'wb') as p:
	pickle.dump(sag_hashes_set, p)
'''
with open(join(save_path, sag_id + '.pkl'), 'rb') as p:
	sag_hashes = pickle.load(p)
	sag_hashes_set = set(sag_hashes)
print('Unpickled SAG L-hashes')
'''
# MG Tetras
mg_id = mg_file.split('/')[-1].split('.')[0]
mg_contigs = get_seqs(mg_file)
# Break up contigs into overlapping subseqs
mg_headers, mg_subs = get_subseqs(mg_contigs, max_contig_len, overlap_len)
mg_tetra_df = pd.DataFrame.from_dict(tetra_cnt(mg_subs))
mg_tetra_df['contig_id'] = mg_headers
mg_tetra_df.set_index('contig_id', inplace=True)
print('MG %s tetranucleotide frequencies calculated' % mg_id)
mg_tetra_df.to_csv(join(save_path, mg_id + '.tsv'), sep='\t')
'''
mg_tetra_df = pd.read_csv(join(save_path, mg_id + '.tsv'), sep='\t', index_col=0,
						header=0)
print('Opened MG tetranucleotide tsv file')
'''
# MG L-seg hash per fragged contig
pass_list = []
for mg_header, mg_frag in zip(mg_headers, mg_subs):
	tmp, mg_Ls = get_subseqs([(mg_header, mg_frag)], 24, 23)
	mg_hashes = calc_seg(mg_Ls)
	mg_hashes.sort(reverse=True)
	mg_hashes_set = set(mg_hashes)
	if sag_hashes_set.intersection(mg_hashes_set):
		print('%s passed identity filter' % mg_header)
		pass_list.append(mg_header)
	else:
		print('%s failed' % mg_header)
with open(join(save_path, mg_id + '.pkl'), 'wb') as p:
	pickle.dump(pass_list, p)
'''
with open(join(save_path, mg_id + '.pkl'), 'rb') as p:
	pass_list = pickle.load(p)
print('Unpickled MG Pass Contigs')
'''
# Set indices for UMAP colorcoding
sag_tetra_df['grouping'] = ['SAG' for x in sag_tetra_df.index]
mg_new_index = []
for index in mg_tetra_df.index:  # TODO: breaking with GCF_001580535.SAG.fasta
	trimmed_index = index.rsplit('_', 1)[0]
	if (index in pass_list) and (trimmed_index in sag_raw_contig_headers):
		mg_new_index.append('TruePos')
	elif (index in pass_list) and (trimmed_index not in sag_raw_contig_headers):
			mg_new_index.append('FalsePos')
	elif (index not in pass_list) and (trimmed_index in sag_raw_contig_headers):
		mg_new_index.append('FalseNeg')
	elif (index not in pass_list) and (trimmed_index not in sag_raw_contig_headers):
		mg_new_index.append('TrueNeg')

mg_tetra_df['grouping'] = mg_new_index

concat_df = pd.concat([sag_tetra_df, mg_tetra_df])
grouping = concat_df['grouping']
sorter = ['TrueNeg', 'TruePos', 'SAG', 'FalseNeg', 'FalsePos']
sorterIndex = dict(zip(sorter,range(len(sorter))))
concat_df['Rank'] = concat_df['grouping'].map(sorterIndex)
print(concat_df.groupby('grouping')[['aaaa']].count())
concat_df.set_index('grouping', inplace=True)
concat_df.sort_values(by=['Rank'], inplace=True)
concat_df.drop(['Rank'], axis=1, inplace=True)

features = concat_df.values
targets = concat_df.index.values
targets_ints = [x[0] for x in enumerate(targets, start=0)]

data = plot_umap(concat_df, save_path, n_neighbors=30, min_dist=0.0,
							n_components=2, random_state=42
							)

gmm = GMM(n_components=5, covariance_type='full', random_state=42).fit(data)
probs = gmm.predict_proba(data)
probs_df = pd.DataFrame(data=probs.round(3), index=grouping)
probs_df.reset_index(inplace=True)
print(probs_df.groupby('grouping').sum())

labels = gmm.predict(data)

plot_gmm(sag_id, save_path, gmm, data, targets)



















'''
pca = PCA(0.99, whiten=True)
data = pca.fit_transform(features)
print(data.shape)
'''
'''
n_components = np.arange(1, 20)
models = [GMM(n, covariance_type='full', random_state=42).fit(data)
		  for n in n_components]
plt.plot(n_components, [m.bic(data) for m in models], label='BIC')
plt.plot(n_components, [m.aic(data) for m in models], label='AIC')
plt.legend(loc='best')
plt.xlabel('n_components');
plot_file_name = sag_id + '.' + 'AIC_BIC.PCA.png'
plot_save_path = join(save_path, plot_file_name)
plt.savefig(plot_save_path, bbox_inches="tight")
plt.clf()
'''




'''
ax = sns.scatterplot(x=data[:, 0], y=data[:, 1], hue=targets)
plt.gca().set_aspect('equal', 'datalim')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plot_file_name = sag_id + '.' + 'GMM_predict.UMAP.png'
plot_save_path = join(save_path, plot_file_name)
plt.savefig(plot_save_path, bbox_inches="tight")
plt.clf()
'''

'''
### Try UMAP ###
plot_umap(concat_df, save_path)
for n in (2, 5, 10, 20, 50, 100, 200):
	plot_umap(concat_df, save_path, n_neighbors=n,
				title='n_neighbors = {}'.format(n)
				)
for d in (0.0, 0.1, 0.25, 0.5, 0.8, 0.99):
	plot_umap(concat_df, save_path, min_dist=d,
				title='min_dist = {}'.format(d)
				)

for m in ('euclidean', 'l2', 'l1', 'manhattan', 'cityblock', 'braycurtis',
			'canberra', 'chebyshev', 'correlation', 'cosine', 'dice',
			'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'matching',
			'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean',
			'sokalmichener', 'sokalsneath', 'sqeuclidean'):
	print(m)
	plot_umap(concat_df, save_path, metric=m,
				title='metric = {}'.format(m)
				)
'''

# Process the SAGs
'''
sag_files = [join(sag_path, f) for f in listdir(sag_path) if isfile(join(sag_path, f))]

sag_tetra_df_list = []
for sag in sag_files:
	sag_id = sag.split('/')[-1].split('.')[0]
	sag_contigs = get_seqs(sag)
	# Break up contigs into overlapping subseqs
	sag_subs = get_subseqs(sag_contigs, max_contig_len, overlap_len)
	sag_tetra_df = pd.DataFrame.from_dict(tetra_cnt(sag_subs))
	sag_tetra_df['sag_id'] = [sag_id for x in sag_tetra_df.index]
	sag_tetra_df.set_index('sag_id', inplace=True)
	sag_tetra_df_list.append(sag_tetra_df)
	print('SAG %s tetranucleotide frequencies calculated' % sag_id)

concat_df = pd.concat(sag_tetra_df_list)
concat_df.to_csv(join(save_path, 'SAG_tetras.tsv'), sep='\t')

concat_df = pd.read_csv(join(save_path, 'SAG_tetras.tsv'), sep='\t', index_col=0,
						header=0)
'''
'''
mapper = umap.UMAP(n_neighbors=136).fit(features, targets_ints)
test_embedding = mapper.transform(features)
classes = list(set(targets))

fig, ax = plt.subplots(1, figsize=(14, 10))
ax = sns.scatterplot(*mapper.embedding_.T, hue=targets)
plt.gca().set_aspect('equal', 'datalim')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.title('Fashion MNIST Train Digits Embedded via UMAP Transform')
plot_file_name = 'Training_data.png'
plot_save_path = join(save_path, plot_file_name)
plt.savefig(plot_save_path, bbox_inches="tight")
plt.clf()


fig, ax = plt.subplots(1, figsize=(14, 10))
ax = sns.scatterplot(*mapper.embedding_.T, hue=targets)
plt.gca().set_aspect('equal', 'datalim')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.title('Fashion MNIST Train Digits Embedded via UMAP')
plot_file_name = 'Test_data.png'
plot_save_path = join(save_path, plot_file_name)
plt.savefig(plot_save_path, bbox_inches="tight")
plt.clf()
'''

# Old code
'''
class KDEClassifier(BaseEstimator, ClassifierMixin):
	"""Bayesian generative classification based on KDE
	
	Parameters
	----------
	bandwidth : float
		the kernel bandwidth within each class
	kernel : str
		the kernel name, passed to KernelDensity
	"""
	def __init__(self, bandwidth=1.0, kernel='gaussian'):
		self.bandwidth = bandwidth
		self.kernel = kernel
		
	def fit(self, X, y):
		self.classes_ = np.sort(np.unique(y))
		training_sets = [X[y == yi] for yi in self.classes_]
		self.models_ = [KernelDensity(bandwidth=self.bandwidth,
									  kernel=self.kernel).fit(Xi)
						for Xi in training_sets]
		self.logpriors_ = [np.log(Xi.shape[0] / X.shape[0])
						   for Xi in training_sets]
		return self
		
	def predict_proba(self, X):
		logprobs = np.array([model.score_samples(X)
							 for model in self.models_]).T
		result = np.exp(logprobs + self.logpriors_)
		return result / result.sum(1, keepdims=True)
		
	def predict(self, X):
		return self.classes_[np.argmax(self.predict_proba(X), 1)]


### KDE ###
bandwidths = 10 ** np.linspace(0, 2, 100)
grid = GridSearchCV(KDEClassifier(), {'bandwidth': bandwidths})
grid.fit(features, targets)

scores = [val.mean_validation_score for val in grid.grid_scores_]
plt.semilogx(bandwidths, scores)
plt.xlabel('bandwidth')
plt.ylabel('accuracy')
plt.title('KDE Model Performance')
plot_save_path = join(save_path, 'KDE_model_perform.png')
plt.savefig(plot_save_path, bbox_inches="tight")
plt.clf()
print(grid.best_params_)
print('accuracy =', grid.cv_results_)
'''
### PCA AIC ###
'''
pca = PCA(0.99, whiten=True)
data = pca.fit_transform(features)
print(data.shape)

n_components = np.arange(20, 150, 10)
models = [GMM(n, covariance_type='full', random_state=0)
		  for n in n_components]
aics = [model.fit(data).aic(data) for model in models]
plt.plot(n_components, aics)
plot_save_path = join(save_path, 'GMM_AIC.png')
plt.savefig(plot_save_path, bbox_inches="tight")
plt.clf()
'''

'''
# Process the MetaG fasta
mg_contigs = []
with open(mg_fasta, 'r') as f:
	data = f.read()
	split_data = data.split('>')
	for reccord in split_data:
		split_rec = reccord.split('\n')
		seq = ''.join(split_rec[1:])
		if seq != '':
			mg_contigs.append(seq)
# Break up contigs into overlapping subseqs
mg_subs = get_subseqs(mg_contigs, max_contig_len, overlap_len)

mg_tetra_df = pd.DataFrame.from_dict(tetra_cnt(mg_subs))
mg_tetra_df['contig_id'] = ['metag0' for x in mg_tetra_df.index]
mg_tetra_df.set_index('contig_id', inplace=True)
print('Metagenome tetranucleotide frequencies calculated')

concat_df = pd.concat([mg_tetra_df, sag0_tetra_df, sag1_tetra_df])
'''
# Don't need to scale tetra Hz data, I think?

# PCA
'''
features = concat_df.values
targets = concat_df.index.values
x = StandardScaler().fit_transform(features)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
pc_df = pd.DataFrame(data = principalComponents
			 , columns = ['principal component 1', 'principal component 2'])
final_df = pc_df.set_index(concat_df.index)
print(final_df.head())

ax = sns.scatterplot(x="principal component 1", y="principal component 2",
						hue=final_df.index, data=final_df)
plt.savefig('standard_pca.png')
plt.clf()
'''