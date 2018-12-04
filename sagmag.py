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
		kmer_list = [''.join(x) for x in get_tetra(clean_seq)]
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
		summed_tetras = 0
		for tetra in dedup_dict.keys():
			if dedup_dict[tetra] != '':
				t_prop = (tmp_dict[tetra] 
							+ tmp_dict[dedup_dict[tetra]]) / total_kmer_cnt
				tetra_prop_dict[tetra] = t_prop
				summed_tetras += t_prop 
			else:
				t_prop = tmp_dict[tetra] / total_kmer_cnt
				tetra_prop_dict[tetra] = t_prop
				summed_tetras += t_prop
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

def get_tetra(seq, n=4):  # found on internet
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result


### Start Main ###
sag_fasta = sys.argv[1]
mg_fasta = sys.argv[2]

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
sag_tetra_df = pd.DataFrame.from_dict(tetra_cnt(sag_contigs))
sag_tetra_df['contig_id'] = ['sag_0' for x in sag_tetra_df.index]
sag_tetra_df.set_index('contig_id', inplace=True)
print(sag_tetra_df.head())
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
mg_tetra_df = pd.DataFrame.from_dict(tetra_cnt(mg_contigs))
mg_tetra_df['contig_id'] = ['contig_0' for x in mg_tetra_df.index]
mg_tetra_df.set_index('contig_id', inplace=True)
print(mg_tetra_df.head())
print('Metagenome tetranucleotide frequencies calculated')

concat_df = pd.concat([sag_tetra_df, mg_tetra_df])
print(concat_df.shape)

# Don't need to scale tetra Hz data, I think?

# PCA
features = concat_df.values
targets = concat_df.index.values
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(features)
pc_df = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2'])
final_df = pc_df.set_index(concat_df.index)
print(final_df.head())

ax = sns.scatterplot(x="principal component 1", y="principal component 2",
						hue=final_df.index, data=final_df)
plt.savefig('output.png')