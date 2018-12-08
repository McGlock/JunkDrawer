import sys
from itertools import islice
import csv
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import pandas as pd
import matplotlib
matplotlib.use('agg')
import seaborn as sns
import matplotlib.pyplot as plt



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
	offset = l_max - o_lap
	for i in range(0, len(seq), offset):
		if i+l_max < len(seq):
			frag = seq[i:i+l_max]
		else:
			frag = seq[-l_max:]
		seq_frags.append(frag)
	#seq_frags = [seq[i:i+l_max] for i in
	#				 range(0, len(seq), l_max)]
	#tail_frag = seq[-l_max:]
	#seq_frags[-1] = tail_frag
	return seq_frags

def get_subseqs(seqs, n, o_lap):
	all_sub_seqs = []
	for seq in seqs:
		clean_seq = seq.strip('\n').lower()
		sub_list = get_frags(clean_seq, n, o_lap)
		all_sub_seqs.extend(sub_list)	
	return all_sub_seqs


def plot_umap(df, n_neighbors=15, min_dist=0.1,
				n_components=2, metric='euclidean', title=''):
	fit = umap.UMAP(
		n_neighbors=n_neighbors,
		min_dist=min_dist,
		n_components=n_components,
		metric=metric
	)
	features = concat_df.values
	targets = concat_df.index.values
	colors = ['r', 'g', 'b']
	color_list = [colors[x[0]] for x in enumerate(set(targets))]

	embedding = fit.fit_transform(features)
	ax = sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], hue=targets)
	plt.gca().set_aspect('equal', 'datalim')
	plt.title(title, fontsize=18)
	plot_file_name = '_'.join([str(n_neighbors), str(min_dist),
								str(n_components), metric]) + '.png'
	plt.savefig(plot_file_name)
	plt.clf()



### Start Main ###
sag_fasta = sys.argv[1]
mg_fasta = sys.argv[2]
max_contig_len = int(sys.argv[3])
overlap_len = int(sys.argv[4])
print('Max contig size is %s bp' % max_contig_len)

# Process the SAG fasta
sag_contigs = []
with open(sag_fasta, 'r') as f:
	data = f.read()
	split_data = data.split('>')
	for reccord in split_data:
		split_rec = reccord.split('\n')
		seq = ''.join(split_rec[1:])
		if seq != '':
			sag_contigs.append(seq)
# Break up contigs into overlapping subseqs
sag_subs = get_subseqs(sag_contigs, max_contig_len, overlap_len)
sag_tetra_df = pd.DataFrame.from_dict(tetra_cnt(sag_subs))
sag_tetra_df['contig_id'] = ['sag_0' for x in sag_tetra_df.index]
sag_tetra_df.set_index('contig_id', inplace=True)
print('SAG tetranucleotide frequencies calculated')
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
mg_tetra_df['contig_id'] = ['contig_0' for x in mg_tetra_df.index]
mg_tetra_df.set_index('contig_id', inplace=True)
print('Metagenome tetranucleotide frequencies calculated')

concat_df = pd.concat([sag_tetra_df, mg_tetra_df])

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

### Try UMAP ###
import numpy as np
from sklearn.datasets import load_iris, load_digits
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sns.set(style='white', context='notebook', rc={'figure.figsize':(14,10)})

import umap

for n in (2, 5, 10, 20, 50, 100, 200):
	plot_umap(concat_df, n_neighbors=n, title='n_neighbors = {}'.format(n))
for d in (0.0, 0.1, 0.25, 0.5, 0.8, 0.99):
    plot_umap(concat_df, min_dist=d, title='min_dist = {}'.format(d))

### NEXT STEPS: Transforming New Data with UMAP ###
### https://umap-learn.readthedocs.io/en/latest/transform.html ###
