import sys
from os import listdir
from os.path import isfile, join, isdir, basename, dirname
from itertools import islice, product, combinations
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import hdbscan
import umap
from sklearn.preprocessing import StandardScaler, MinMaxScaler, normalize
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
import functools


sns.set(style='white', context='notebook', rc={'figure.figsize':(14,10)})


nuc_naughty_list = ['r', 'y', 's', 'w', 'k',
					'm', 'b', 'd', 'h', 'v', 'n']


def tetra_cnt(seq_list):
	# Dict of all tetramers
	tetra_cnt_dict = {''.join(x):[] for x in product('atgc', repeat=4)}
	# count up all tetramers and also populate the tetra dict
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


def get_kmer(seq, n):
	"Returns a sliding window (of width n) over data from the iterable"
	"   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...				   "
	it = iter(seq)
	result = tuple(islice(it, n))
	if len(result) == n:
		yield result
	for elem in it:
		result = result[1:] + (elem,)
		yield result


def get_frags(seq, l_max, o_lap):
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
			#if any(nuc in frag for nuc in nuc_naughty_list):
			#	1+1  # TODO: not sure how to handle these
			#else:
			seq_frags.append(frag)
	#elif 'n' not in seq:
	else:
		seq_frags.append(seq)
	#else:
	#	1+1  # TODO: not sure how to handle these
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
			if nuc not in nuc_naughty_list:  # TODO: don't know what to do with these
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
		nsig_width = nsig * width
		nsig_height = nsig * height
		ax.add_patch(Ellipse(position, nsig_width, nsig_height,
							 angle, **kwargs))
	return position, nsig_width, nsig_height, angle

	

def plot_gmm(sg_id, sv_pth, gmm, X, targets, label=True, ax=None):
	# Create scatter with error group labels
	ax = sns.scatterplot(x=X[:, 0], y=X[:, 1], hue=targets)
	plt.gca().set_aspect('equal', 'datalim')
	w_factor = 0.1 / gmm.weights_.max()
	for pos, covar, w in zip(gmm.means_, gmm.covariances_, gmm.weights_):
		draw_ellipse(pos, covar, alpha=w * w_factor)
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plot_file_name = sag_id + '.' + 'GMM_error_classes.UMAP.png'
	plot_save_path = join(save_path, plot_file_name)
	plt.savefig(plot_save_path, bbox_inches="tight")
	plt.clf()
	# Predict labels based on GMM
	labels = gmm.predict(X)
	# Get associated probabilities for predicted labels
	probs = gmm.predict_proba(X)
	probs_df = pd.DataFrame(data=probs.round(3), index=targets)
	probs_df.reset_index(inplace=True)
	# Create scatter with predicted membership labels
	if label:
		ax = sns.scatterplot(x=X[:, 0], y=X[:, 1], hue=labels)
		plt.gca().set_aspect('equal', 'datalim')
		w_factor = 0.1 / gmm.weights_.max()
		for pos, covar, w in zip(gmm.means_, gmm.covariances_, gmm.weights_):
			draw_ellipse(pos, covar, alpha=w * w_factor)
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plot_file_name = sag_id + '.' + 'GMM_predicted_labels.UMAP.png'
		plot_save_path = join(save_path, plot_file_name)
		plt.savefig(plot_save_path, bbox_inches="tight")
		plt.clf()
		probs_df['p_labels'] = labels

	return probs_df


def plot_ellispe_membership(df, plot_save_path, mean, covar):
	# Draw ellispe that colors by membership
	position, width, height, angle = draw_ellipse(mean, covar,
													edgecolor='green',
													linewidth=1, alpha=0.1
													)
	cos_angle = np.cos(np.radians(180.-angle))
	sin_angle = np.sin(np.radians(180.-angle))
	xc = df[df.columns[0]] - position[0]
	yc = df[df.columns[1]] - position[1]
	xct = xc * cos_angle - yc * sin_angle
	yct = xc * sin_angle + yc * cos_angle 
	rad_cc = (xct**2/(width/2.)**2) + (yct**2/(height/2.)**2)
	membr_category = []
	for r in rad_cc:
		if r <= 1.:
			# point in ellipse
			membr_category.append('SAG+')
		else:
			# point not in ellipse
			membr_category.append('Not SAG')
	if membr_category[0] == 'SAG+':
		pal = ['green', 'gray']
	else:
		pal = ['gray', 'green']
	ax = sns.scatterplot(x=df[df.columns[0]], y=df[df.columns[1]],
							hue=membr_category,	palette=pal
							)
	plt.gca().set_aspect('equal', 'datalim')
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.savefig(plot_save_path, bbox_inches="tight")
	plt.clf()
	# add membership to df
	isSAG_col = '_'.join(['isSAG', df.columns[0], df.columns[1]])
	df[isSAG_col] = [1 if x == 'SAG+' else 0 for x in membr_category]
	df.drop(df.columns[:2], axis=1, inplace=True)

	return df, isSAG_col


def plot_ellispe_error(df, plot_save_path, mean, covar):
	# Draw ellispe that colors by error stats
	ax = sns.scatterplot(x=df[df.columns[0]], y=df[df.columns[1]], hue=df.index)
	draw_ellipse(mean, covar, alpha=0.1)
	plt.gca().set_aspect('equal', 'datalim')
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.savefig(plot_save_path, bbox_inches="tight")
	plt.clf()


def calc_err(df):
	# build error type df
	idf_cnt_df = df.groupby('idf_errors')[df.columns[0]].count().reset_index()
	idf_cnt_df.columns = ['err_type', 'idf_errors']
	idf_TP = idf_cnt_df.loc[idf_cnt_df['err_type'] == 'TruePos', 'idf_errors'].values[0]
	idf_FP = idf_cnt_df.loc[idf_cnt_df['err_type'] == 'FalsePos', 'idf_errors'].values[0]
	idf_FN = idf_cnt_df.loc[idf_cnt_df['err_type'] == 'FalseNeg', 'idf_errors'].values[0]
	idf_TN = idf_cnt_df.loc[idf_cnt_df['err_type'] == 'TrueNeg', 'idf_errors'].values[0]
	
	idf_cnt_df['idf_precision'] = idf_TP/(idf_TP + idf_FP)
	idf_cnt_df['idf_sensitivity'] = idf_TP/(idf_TP + idf_FN)
	idf_cnt_df['idf_specificity'] = idf_TN/(idf_TN + idf_FP)

	thf_cnt_df = df.groupby('thf_errors')[df.columns[0]].count().reset_index()
	thf_cnt_df.columns = ['err_type', 'thf_errors']
	thf_TP = thf_cnt_df.loc[thf_cnt_df['err_type'] == 'TruePos', 'thf_errors'].values[0]
	thf_FP = thf_cnt_df.loc[thf_cnt_df['err_type'] == 'FalsePos', 'thf_errors'].values[0]
	thf_FN = thf_cnt_df.loc[thf_cnt_df['err_type'] == 'FalseNeg', 'thf_errors'].values[0]
	thf_TN = thf_cnt_df.loc[thf_cnt_df['err_type'] == 'TrueNeg', 'thf_errors'].values[0]
	
	thf_cnt_df['thf_precision'] = thf_TP/(thf_TP + thf_FP)
	thf_cnt_df['thf_sensitivity'] = thf_TP/(thf_TP + thf_FN)
	thf_cnt_df['thf_specificity'] = thf_TN/(thf_TN + thf_FP)

	htf_cnt_df = df.groupby('htf_errors')[df.columns[0]].count().reset_index()
	htf_cnt_df.columns = ['err_type', 'htf_errors']
	htf_TP = htf_cnt_df.loc[htf_cnt_df['err_type'] == 'TruePos', 'htf_errors'].values[0]
	htf_FP = htf_cnt_df.loc[htf_cnt_df['err_type'] == 'FalsePos', 'htf_errors'].values[0]
	htf_FN = htf_cnt_df.loc[htf_cnt_df['err_type'] == 'FalseNeg', 'htf_errors'].values[0]
	htf_TN = htf_cnt_df.loc[htf_cnt_df['err_type'] == 'TrueNeg', 'htf_errors'].values[0]
	
	htf_cnt_df['htf_precision'] = htf_TP/(htf_TP + htf_FP)
	htf_cnt_df['htf_sensitivity'] = htf_TP/(htf_TP + htf_FN)
	htf_cnt_df['htf_specificity'] = htf_TN/(htf_TN + htf_FP)

	err_df_list = [idf_cnt_df, thf_cnt_df, htf_cnt_df]

	#error_df = pd.merge(idf_cnt_df, thf_cnt_df, on='err_type')
	error_df = functools.reduce(lambda left,right: pd.merge(left,right,on='err_type'),
									err_df_list
									)

	return error_df


def mock_SAG(fasta_file):
	# currently just returns half of the genome as a mock SAG
	genome_contigs = get_seqs(fasta_file)
	if len(genome_contigs) != 1:
		half_list = genome_contigs[::2]
	else:
		header = genome_contigs[0][0]
		seq = genome_contigs[0][1]
		half_list = [(header,seq[:int(len(seq)/2)])]
	all_headers = [x[0] for x in genome_contigs]
	return half_list, all_headers


def main():
	sag_path = sys.argv[1]
	mg_file = sys.argv[2]
	sag_abund_file = sys.argv[3]
	mg_abund_file = sys.argv[4]
	max_contig_len = int(sys.argv[5])
	overlap_len = int(sys.argv[6])
	save_path = sys.argv[7]
	num_components = int(sys.argv[8])

	# magic numbers
	#num_components = 2

	# TODO: build argv interface
	# TODO: extract all magic numbers to be built into argv
	# TODO: Refactor all definitions, they are a mess :/
	# TODO: Fit UMAP output to model and test with AIC, which model fits best?
	# TODO: combine error into one value, remove subseq map file creation
	# TODO: add abundance(coverage) data to analysis

	if isdir(sag_path):
		print('[SAG+]: Directory specified, looking for .fasta files within')
		sag_list = [join(sag_path, f) for f in
								listdir(sag_path) if '.fasta' in f
								]
	elif isfile(sag_path):
		print('[SAG+]: File specified, processing %s' % basename(sag_path))
		sag_list = [sag_path]

	error_df_list = []
	subseq_map_list = []
	for sag_file in sag_list:
		sag_basename = basename(sag_file)
		sag_id = sag_basename.split('.')[0]
		
		# SAG Tetras
		if isfile(join(save_path, sag_id + '.tsv')):
			#sag_contigs, sag_raw_contig_headers = mock_SAG(sag_file)
			sag_tetra_df = pd.read_csv(join(save_path, sag_id + '.tsv'),
										sep='\t', index_col=0, header=0)

			with open(join(save_path, sag_id + '.headers.pkl'), 'rb') as p:
				sag_raw_contig_headers = pickle.load(p)
			print('[SAG+]: Found %s SAG tetranucleotide tsv file' % sag_id)
		else:
			print('[SAG+]: Calculating tetramer frequencies for %s' % sag_id)
			### Used for Mock SAGs (need to change when done testing)
			sag_contigs, sag_raw_contig_headers = mock_SAG(sag_file)
			#sag_contigs = get_seqs(sag_file)
			sag_headers, sag_subs = get_subseqs(sag_contigs, max_contig_len, overlap_len)
			sag_tetra_df = pd.DataFrame.from_dict(tetra_cnt(sag_subs))
			sag_tetra_df['contig_id'] = sag_headers
			sag_tetra_df.set_index('contig_id', inplace=True)
			sag_tetra_df.to_csv(join(save_path, sag_id + '.tsv'), sep='\t')
			with open(join(save_path, sag_id + '.headers.pkl'), 'wb') as p:
				pickle.dump(sag_raw_contig_headers, p)

		'''
		# SAG coverage info
		sag_abund_df = pd.read_csv(sag_abund_file, header=0, sep='\t')
		sag_tot_depth_dict = {x[0]: sag_abund_df.loc[sag_abund_df['contigName']
								== x[0]]['totalAvgDepth'].values[0]
								for x in sag_contigs
								}
		'''

		# SAG subseqs kmer hashing
		if isfile(join(save_path, sag_id + '.pkl')):
			with open(join(save_path, sag_id + '.pkl'), 'rb') as p:
				sag_hashes = pickle.load(p)
				sag_hashes_set = set(sag_hashes)
			print('[SAG+]: Unpickled %s kmer hashes' % sag_id)
		else:
			print('[SAG+]: Calculating kmer hashes for %s' % sag_id)
			tmp, sag_Ls = get_subseqs(sag_contigs, 24, 23)
			sag_hashes = calc_seg(sag_Ls)
			sag_hashes.sort(reverse=True)
			sag_hashes_set = set(sag_hashes)
			with open(join(save_path, sag_id + '.pkl'), 'wb') as p:
				pickle.dump(sag_hashes_set, p)
		
		# MG Tetras
		mg_basename = basename(mg_file)
		mg_id = mg_basename.split('.')[0]
		if isfile(join(save_path, mg_id + '.tsv')):
			mg_tetra_df = pd.read_csv(join(save_path, mg_id + '.tsv'), sep='\t', index_col=0,
									header=0)
			mg_contigs = get_seqs(mg_file)
			#mg_contigs_trm_header = [(x[0].rsplit('|', 1)[0], x[1]) for x in mg_contigs]
			mg_headers, mg_subs = get_subseqs(mg_contigs, max_contig_len, overlap_len)
			print('[SAG+]: Found %s MetaG tetranucleotide tsv file' % mg_id)
		else:
			print('[SAG+]: Calculating tetramer frequencies for %s' % mg_id)
			mg_contigs = get_seqs(mg_file)
			#mg_contigs_trm_header = [(x[0].rsplit('|', 1)[1], x[1]) for x in mg_contigs]
			mg_headers, mg_subs = get_subseqs(mg_contigs, max_contig_len, overlap_len)
			mg_tetra_df = pd.DataFrame.from_dict(tetra_cnt(mg_subs))
			mg_tetra_df['contig_id'] = mg_headers
			mg_tetra_df.set_index('contig_id', inplace=True)
			mg_tetra_df.to_csv(join(save_path, mg_id + '.tsv'), sep='\t')

		'''
		# MetaG coverage info
		mg_abund_df = pd.read_csv(mg_abund_file, header=0, sep='\t')
		mg_tot_depth_dict = {x[0].rsplit('|', 1)[0]: mg_abund_df.loc[mg_abund_df['contigName']
								== x[0].rsplit('|', 1)[0]]['totalAvgDepth'].values[0]
								for x in mg_contigs
								}
		'''

		# MG subseqs L-mer hash, compare to SAG hashes
		if isfile(join(save_path, sag_id + '.kmer_recruit.pkl')): 
			with open(join(save_path, sag_id + '.kmer_recruit.pkl'), 'rb') as p:
				pass_list = pickle.load(p)
			print('[SAG+]: Unpickled %s kmer ID filter' % sag_id)
		else:
			print('[SAG+]: Performing kmer ID filtering')
			pass_list = []
			for mg_header, mg_frag in zip(mg_headers, mg_subs):  # TODO: this is really slow :(
				tmp, mg_Ls = get_subseqs([(mg_header, mg_frag)], 24, 23)
				mg_hashes = calc_seg(mg_Ls)
				mg_hashes.sort(reverse=True)
				mg_hashes_set = set(mg_hashes)
				if sag_hashes_set.intersection(mg_hashes_set):
					#print('%s passed identity filter' % mg_header)
					pass_list.append(mg_header)
				#else:
					#print('%s failed' % mg_header)
			with open(join(save_path, sag_id + '.kmer_recruit.pkl'), 'wb') as p:
				pickle.dump(pass_list, p)
		
		### Used for seq tracking and error analysis
		# Look at ID filter (idf) error types
		mg_idf_errors = []
		for index in mg_tetra_df.index:
			trimmed_index = index.rsplit('|', 1)[1].rsplit('_', 1)[0]
			if (index in pass_list) and (trimmed_index in sag_raw_contig_headers):
				mg_idf_errors.append('TruePos')
			elif (index in pass_list) and (trimmed_index not in sag_raw_contig_headers):
					mg_idf_errors.append('FalsePos')
			elif (index not in pass_list) and (trimmed_index in sag_raw_contig_headers):
				mg_idf_errors.append('FalseNeg')
			elif (index not in pass_list) and (trimmed_index not in sag_raw_contig_headers):
				mg_idf_errors.append('TrueNeg')
		track_recruits_df = pd.DataFrame(mg_idf_errors, columns=['idf_errors'],
											index=mg_tetra_df.index)
		# Set MetaG contig index to genome id for error tracking
		#mg_contig_index = [x.rsplit('|', 1)[1] for x in mg_tetra_df.index]
		#mg_tetra_df.reset_index(inplace=True)
		#mg_tetra_df['contig_id'] = mg_contig_index
		#mg_tetra_df.set_index('contig_id', inplace=True)

		concat_tetra_df = pd.concat([sag_tetra_df, mg_tetra_df])
		normed_tetra_df = pd.DataFrame(normalize(concat_tetra_df.values),
										columns=concat_tetra_df.columns,
										index=concat_tetra_df.index
										)
		sag_tetra_df['data_type'] = ['SAG' for x in sag_tetra_df.index]
		mg_tetra_df['data_type'] = ['MG' for x in mg_tetra_df.index]

		#sorter = ['TrueNeg', 'TruePos', 'SAG', 'FalseNeg', 'FalsePos']
		#sorterIndex = dict(zip(sorter,range(len(sorter))))
		#concat_df['Rank'] = concat_df['idf_errors'].map(sorterIndex)
		#concat_df.sort_values(by=['Rank'], inplace=True)
		#sorted_subseq_ids = concat_df.index.values
		#idf_df = concat_df.set_index('idf_errors')
		#idf_df.drop(['Rank'], axis=1, inplace=True)
		#idf_df = pd.DataFrame(normalize(idf_df.values), columns=idf_df.columns,
		#						index=idf_df.index)
		### END

		features = normed_tetra_df.values
		targets = normed_tetra_df.index.values
		targets_ints = [x[0] for x in enumerate(targets, start=0)]

		print('[SAG+]: Dimension reduction with UMAP')
		#data = plot_umap(idf_df, save_path, n_neighbors=30, min_dist=0.0,
		#							n_components=num_components, random_state=42,
		#							metric='manhattan'
		#							)
		
		fit = umap.UMAP(n_neighbors=30, min_dist=0.0,
						n_components=num_components, metric='manhattan',
						random_state=42
						)
		data = fit.fit_transform(features)
		pc_col_names = ['pc' + str(x) for x in range(1, num_components + 1)]
		umap_df = pd.DataFrame(data, columns=pc_col_names, index=targets)
		# build SAG GMM for H0 test file (htf)
		print('[SAG+]: Calculating SAG GMMs')
		sag_umap_df = umap_df.loc[umap_df.index.isin(sag_tetra_df.index)]
		mg_umap_df = umap_df.loc[umap_df.index.isin(mg_tetra_df.index)]
		n_components = np.arange(1, 100, 1)
		models = [GMM(n, covariance_type='tied', random_state=42)
			  for n in n_components]
		bics = [model.fit(sag_umap_df.values).bic(sag_umap_df.values) for model in models]
		aics = [model.fit(sag_umap_df.values).aic(sag_umap_df.values) for model in models]
		min_bic_comp = n_components[bics.index(min(bics))]
		min_aic_comp = n_components[aics.index(min(aics))]
		print(min_bic_comp, min_aic_comp)
		gmm = GMM(n_components=min_bic_comp, covariance_type='tied',
					random_state=42).fit(sag_umap_df.values)
		sag_scores = gmm.score_samples(sag_umap_df.values)
		sag_scores_df = pd.DataFrame(data=sag_scores, index=sag_umap_df.index)
		sag_score_min = min(sag_scores_df.values)[0]
		sag_score_max = max(sag_scores_df.values)[0]

		mg_scores = gmm.score_samples(mg_umap_df.values)
		mg_scores_df = pd.DataFrame(data=mg_scores, index=mg_umap_df.index)
		mg_score_min = min(mg_scores_df.values)[0]
		mg_score_max = max(mg_scores_df.values)[0]
		mg_htf_errors = []
		for index, row in mg_scores_df.iterrows():
			score = row[0]
			trimmed_header = index.rsplit('|', 1)[1].rsplit('_', 1)[0]
			if ((sag_score_min <= score <= sag_score_max) and
					(trimmed_header in sag_raw_contig_headers)):
				mg_htf_errors.append('TruePos')
			elif ((sag_score_min <= score <= sag_score_max) and
					(trimmed_header not in sag_raw_contig_headers)):
				mg_htf_errors.append('FalsePos')
			elif (((sag_score_min > score) or (score > sag_score_max)) and
					(trimmed_header in sag_raw_contig_headers)):
				mg_htf_errors.append('FalseNeg')
			elif (((sag_score_min > score) or (score > sag_score_max)) and
					(trimmed_header not in sag_raw_contig_headers)):
				mg_htf_errors.append('TrueNeg')

		pc_pair_error_df_list = []
		pc_pair_subseq_map_list = []
		isSAG_cols = []
		# TODO: refactoring this loop
		for pc_pair in combinations(pc_col_names, 2):
			pc_pair = list(pc_pair)
			subset_df = umap_df[pc_pair]
			subsag_umap_df = sag_umap_df[pc_pair]
			submg_umap_df = mg_umap_df[pc_pair]

			sag_std = subsag_umap_df.std().values
			sag_mean = subsag_umap_df.mean().values
			#sag_covar = subsag_umap_df.cov().values
			sag_corr = subsag_umap_df.corr().values
			### Used for seq tracking and error analysis
			# Draw ellispe that colors by L-mer error stats
			print('[SAG+]: Plotting clusting with kmer error stats')
			print('[SAG+]: Including SAG correlation distribution ellispe')
			error_file_name = '.'.join([sag_id, '_'.join(pc_pair),
										'UMAP_Error_Ellipse', 'png'])
			error_save_path = join(save_path, error_file_name)
			plot_ellispe_error(submg_umap_df, error_save_path, sag_mean, sag_corr)
			### END

			# Draw ellispe that colors by membership
			print('[SAG+]: Plotting clusting with subcontig membership')
			print('[SAG+]: Including SAG correlation distribution ellispe')
			memb_file_name = '.'.join([sag_id, '_'.join(pc_pair),
										'UMAP_Ellipse_membership', 'png'])
			memb_save_path = join(save_path, memb_file_name)
			membership_df, isSAG_col = plot_ellispe_membership(submg_umap_df, memb_save_path,
													sag_mean, sag_corr)
			isSAG_cols.append(isSAG_col)
			# add subseq mapping
			membership_df['subseq_header'] = mg_umap_df.index
			
			### Used for seq tracking and error analysis
			# preserve idf errors
			membership_df['idf_errors'] = membership_df.index
			# Look at tetramer Hz filter (thf) error types
			print('[SAG+]: Building error type dataframe')
			mg_thf_errors = []
			for index, row in membership_df.iterrows():
				header = row['subseq_header']
				isSAG = row[isSAG_col]
				trimmed_header = header.rsplit('_', 1)[0]
				if (isSAG == 1) and (trimmed_header in sag_raw_contig_headers):
					mg_thf_errors.append('TruePos')
				elif (isSAG == 1) and (trimmed_header not in sag_raw_contig_headers):
						mg_thf_errors.append('FalsePos')
				elif (isSAG != 1) and (trimmed_header in sag_raw_contig_headers):
					mg_thf_errors.append('FalseNeg')
				elif (isSAG != 1) and (trimmed_header not in sag_raw_contig_headers):
					mg_thf_errors.append('TrueNeg')
			membership_df['thf_errors'] = mg_thf_errors
			membership_df['htf_errors'] = mg_htf_errors
			error_df = calc_err(membership_df)

			membership_df.drop(columns=['idf_errors', 'thf_errors', 'htf_errors'],
								axis=1, inplace=True)
			membership_df['dimension'] = '_'.join(pc_pair)
			error_df['sag_id'] = sag_id
			error_df.set_index('sag_id', inplace=True)
			error_df['dimension'] = '_'.join(pc_pair)
			pc_pair_error_df_list.append(error_df)
			### END
			pc_pair_subseq_map_list.append(membership_df)
		pc_pair_subseq_df = functools.reduce(lambda x, y: pd.merge(  # TODO: Eating all the RAMs, find a fix
															x, y, on=['subseq_header']),
															pc_pair_subseq_map_list
															)
		pc_pair_subseq_df.set_index('subseq_header', inplace=True)
		pc_pair_subseq_df.to_csv(join(save_path, sag_id + '_subseq_map.tsv'), sep='\t')

		### Used for seq tracking and error analysis
		pc_pair_err_df = pd.concat(pc_pair_error_df_list)
		pc_pair_err_df.to_csv(join(save_path, sag_id + '_error_stats.tsv'), sep='\t')
		### END

		error_df_list.append(pc_pair_err_df)
		subseq_map_list.append(pc_pair_subseq_df)

		# get all predicted SAGs
		#SAG_pred_dict = {}
		SAG_pred_list = []
		for index, row in pc_pair_subseq_df.iterrows():
			isSAG_val_list = [row[v] for v in isSAG_cols]
			if sum(isSAG_val_list) == len(isSAG_cols):
				SAG_pred_list.append(index)
		print('[SAG+]: Predicted %s subcontigs for SAG %s' % (str(len(set(SAG_pred_list))), 
																sag_id)
																)

		'''
			seq_header = index.rsplit('_', 1)[0]
			if seq_header in SAG_pred_dict.keys():
				SAG_pred_dict[seq_header].append(sum(isSAG_val_list))
			else:
				SAG_pred_dict[seq_header] = [sum(isSAG_val_list)]
		'''

		# Save Predicted SAG contigs to a fasta file
		fasta_out_file = join(save_path, sag_id + '.predicted_contigs.fasta')
		with open(fasta_out_file, 'w') as fasta_out:
			for header, seq in zip(mg_headers, mg_subs):
				if header in SAG_pred_list:
					fasta_out.write('\n'.join([header, seq]) + '\n')
		print('[SAG+]: Predicted subcontigs saved to %s' % basename(fasta_out_file))

		'''
		# Get only contigs where all subcontigs are in the cluster
		for key in SAG_pred_dict.keys():
			sum_all = sum(SAG_pred_dict[key])
			if sum_all/len(SAG_pred_dict[key]) == len(isSAG_cols):
				1+1#print(key)
			elif sum_all/len(SAG_pred_dict[key]) != 0:
				print(sum_all/len(SAG_pred_dict[key]))
				print(key)
				print(SAG_pred_dict[key])
		'''
		print('[SAG+]: Completed analysis of %s and %s' % (sag_id, mg_id))
	final_subseq_df = pd.concat(subseq_map_list)
	final_subseq_df.to_csv(join(save_path, 'total_subseq_map.tsv'), sep='\t')
	### Used for seq tracking and error analysis
	final_err_df = pd.concat(error_df_list)  # TODO: combine error stats into one value
	final_err_df.to_csv(join(save_path, 'total_error_stats.tsv'), sep='\t')
	### END


if __name__ == "__main__":
	main()




################################################################
### ANYTHING BELOW THIS IS JUNK, WILL BE NUKED AT SOME POINT ###
################################################################



# Old code

'''
gmm = GMM(n_components=5, covariance_type='full', random_state=42).fit(data, targets)
probs = gmm.predict_proba(data)
probs_df = pd.DataFrame(data=probs.round(3), index=grouping)
probs_df.reset_index(inplace=True)
print(probs_df.groupby('grouping').sum())

labels = gmm.predict(data)

plot_gmm(sag_id, save_path, gmm, data, targets)

'''
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
grid = GridSearchCV(KernelDensity(kernel='gaussian'),
					{'bandwidth': bandwidths}
					)
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

'''

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