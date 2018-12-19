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


def plot_ellispe_error(df, error_hue, plot_save_path, mean, covar):
	# Draw ellispe that colors by error stats
	ax = sns.scatterplot(x=df[df.columns[0]], y=df[df.columns[1]], hue=error_hue)
	draw_ellipse(mean, covar, alpha=0.1)
	plt.gca().set_aspect('equal', 'datalim')
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.savefig(plot_save_path, bbox_inches="tight")
	plt.clf()


def calc_err(df):
	# build error type df
	err_df_list = []
	for col in df.columns:
		val_cnt = df[col].value_counts()
		cnt_df = val_cnt.rename_axis('err_type').to_frame(col).reset_index()

		if 'TruePos' in list(cnt_df['err_type']):
			TP = cnt_df.loc[cnt_df['err_type'] == 'TruePos', col].values[0]
		else:
			TP = 0
		if 'FalsePos' in list(cnt_df['err_type']):
			FP = cnt_df.loc[cnt_df['err_type'] == 'FalsePos', col].values[0]
		else:
			FP = 0
		if 'FalseNeg' in list(cnt_df['err_type']):
			FN = cnt_df.loc[cnt_df['err_type'] == 'FalseNeg', col].values[0]
		else:
			FN = 0
		if 'TrueNeg' in list(cnt_df['err_type']):
			TN = cnt_df.loc[cnt_df['err_type'] == 'TrueNeg', col].values[0]
		else:
			TN = 0

		cnt_df[col + '_precision'] = TP/(TP + FP)
		cnt_df[col + '_sensitivity'] = TP/(TP + FN)
		cnt_df[col + '_specificity'] = TN/(TN + FP)

		err_df_list.append(cnt_df)
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
	tax_map = sys.argv[8]
	taxid_file = sys.argv[9]
	num_components = 3 # int(sys.argv[8])

	# magic numbers
	#num_components = 2

	# TODO: build argv interface
	# TODO: extract all magic numbers to be built into argv
	# TODO: Refactor all definitions, they are a mess :/
	# TODO: Fit UMAP output to model and test with AIC, which model fits best?
	# TODO: combine error into one value, remove subseq map file creation
	# TODO: add abundance(coverage) data to analysis
	# TODO: output number of tetras/kmers skipped b/c of ambiguous nucs

	if isdir(sag_path):
		print('[SAG+]: Directory specified, looking for .fasta files within')
		sag_list = [join(sag_path, f) for f in
								listdir(sag_path) if f.split('.')[-1] == 'fasta'
								]
	elif isfile(sag_path):
		print('[SAG+]: File specified, processing %s' % basename(sag_path))
		sag_list = [sag_path]

	error_df_list = []
	taxa_tracking_list = []
	for sag_file in sag_list:
		sag_basename = basename(sag_file)
		sag_id = sag_basename.split('.')[0]
		
		# SAG Tetras
		if isfile(join(save_path, sag_id + '.tsv')):
			sag_contigs, sag_raw_contig_headers = mock_SAG(sag_file)
			sag_tetra_df = pd.read_csv(join(save_path, sag_id + '.tsv'),
										sep='\t', index_col=0, header=0)

			with open(join(save_path, sag_id + '.headers.pkl'), 'rb') as p:
				sag_raw_contig_headers = pickle.load(p)
			print('[SAG+]: Found %s SAG tetranucleotide tsv file' % sag_id)
		else:
			print('[SAG+]: Calculating tetramer frequencies for %s' % sag_id)
			### Used for Mock SAGs (need to change when done testing)
			sag_contigs, sag_raw_contig_headers = mock_SAG(sag_file)
			sag_headers, sag_subs = get_subseqs(sag_contigs, max_contig_len, overlap_len)
			sag_tetra_df = pd.DataFrame.from_dict(tetra_cnt(sag_subs))
			sag_tetra_df['contig_id'] = sag_headers
			sag_tetra_df.set_index('contig_id', inplace=True)
			sag_tetra_df.to_csv(join(save_path, sag_id + '.tsv'), sep='\t')
			with open(join(save_path, sag_id + '.headers.pkl'), 'wb') as p:
				pickle.dump(sag_raw_contig_headers, p)

		# Save Mock SAG for error analysis
		sag_mock_out = join(save_path, sag_id + '.mSAG.fasta')
		with open(sag_mock_out, 'w') as s_out:
			for header, seq in sag_contigs:
				new_header = '>' + header
				s_out.write('\n'.join([new_header, seq]) + '\n')

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
					pass_list.append(mg_header)
			with open(join(save_path, sag_id + '.kmer_recruit.pkl'), 'wb') as p:
				pickle.dump(pass_list, p)
		
		error_dict = {}

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
		
		error_dict['idf_errors'] = mg_idf_errors

		concat_tetra_df = pd.concat([sag_tetra_df, mg_tetra_df])
		normed_tetra_df = pd.DataFrame(normalize(concat_tetra_df.values),
										columns=concat_tetra_df.columns,
										index=concat_tetra_df.index
										)
		sag_tetra_df['data_type'] = ['SAG' for x in sag_tetra_df.index]
		mg_tetra_df['data_type'] = ['MG' for x in mg_tetra_df.index]
		### END

		features = normed_tetra_df.values
		targets = normed_tetra_df.index.values
		targets_ints = [x[0] for x in enumerate(targets, start=0)]

		print('[SAG+]: Dimension reduction with UMAP')
		fit = umap.UMAP(n_neighbors=30, min_dist=0.0,
						n_components=num_components, metric='manhattan',
						random_state=42
						)
		data = fit.fit_transform(features)
		pc_col_names = ['pc' + str(x) for x in range(1, num_components + 1)]
		umap_df = pd.DataFrame(data, columns=pc_col_names, index=targets)
		# build SAG GMM for H0 test file (htf)
		print('[SAG+]: Calculating AIC/BIC')
		sag_umap_df = umap_df.loc[umap_df.index.isin(sag_tetra_df.index)]
		mg_umap_df = umap_df.loc[umap_df.index.isin(mg_tetra_df.index)]
		n_components = np.arange(1, 100, 1)
		models = [GMM(n, covariance_type='tied', random_state=42)
			  for n in n_components]
		bics = [model.fit(sag_umap_df.values).bic(sag_umap_df.values) for model in models]
		aics = [model.fit(sag_umap_df.values).aic(sag_umap_df.values) for model in models]
		min_bic_comp = n_components[bics.index(min(bics))]
		min_aic_comp = n_components[aics.index(min(aics))]
		print('[SAG+]: Min AIC/BIC at %s/%s, respectively' % (min_aic_comp, min_bic_comp))
		print('[SAG+]: Using AIC as guide for GMM components')
		print('[SAG+]: Training GMM on SAG tetras')
		gmm = GMM(n_components=min_aic_comp, covariance_type='tied',
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

		error_dict['htf_errors'] = mg_htf_errors

		pc_pair_error_df_list = []
		pc_pair_subseq_map_list = []
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
			plot_ellispe_error(submg_umap_df, mg_idf_errors, error_save_path,
								sag_mean, sag_corr
								)
			### END

			# Draw ellispe that colors by membership
			print('[SAG+]: Plotting clusting with subcontig membership')
			print('[SAG+]: Including SAG correlation distribution ellispe')
			memb_file_name = '.'.join([sag_id, '_'.join(pc_pair),
										'UMAP_Ellipse_membership', 'png'])
			memb_save_path = join(save_path, memb_file_name)
			membership_df, isSAG_col = plot_ellispe_membership(submg_umap_df,
																memb_save_path,
																sag_mean, sag_corr
																)
			### Used for seq tracking and error analysis
			# Look at tetramer Hz filter (thf) error types
			mg_thf_errors = []
			for index, row in membership_df.iterrows():
				isSAG = row[isSAG_col]
				trimmed_header = index.rsplit('|', 1)[1].rsplit('_', 1)[0]
				if (isSAG == 1) and (trimmed_header in sag_raw_contig_headers):
					mg_thf_errors.append('TruePos')
				elif (isSAG == 1) and (trimmed_header not in sag_raw_contig_headers):
						mg_thf_errors.append('FalsePos')
				elif (isSAG != 1) and (trimmed_header in sag_raw_contig_headers):
					mg_thf_errors.append('FalseNeg')
				elif (isSAG != 1) and (trimmed_header not in sag_raw_contig_headers):
					mg_thf_errors.append('TrueNeg')
			error_dict['_'.join(['thf', '_'.join(pc_pair)])] = mg_thf_errors

		print('[SAG+]: Building error type dataframe')
		all_isSAG_df = pd.DataFrame(error_dict, index=mg_umap_df.index)
		error_df = calc_err(all_isSAG_df)
		error_df['sag_id'] = sag_id
		error_df.set_index('sag_id', inplace=True)
		error_df.to_csv(join(save_path, sag_id + '_error_stats.tsv'), sep='\t')
		### END
		error_df_list.append(error_df)

		# get subcontigs predicted as SAG+ with all filters
		SAG_pred_list = []
		#TP_list = []
		for index, row in all_isSAG_df.iterrows():
			isSAG_sum = sum([1 for x in row if 'Pos' in x])
			#TP_sum = sum([1 for x in row if 'TruePos' in x])
			#if TP_sum == len(all_isSAG_df.columns):
			#	TP_list.append(index)
			if isSAG_sum == len(all_isSAG_df.columns):
				SAG_pred_list.append(index)
		print('[SAG+]: Predicted %s subcontigs for SAG %s' % (str(len(set(SAG_pred_list))), 
																sag_id)
																)
		#print('[SAG+]: Of those, %s were True Positives' % (str(len(set(TP_list)))))
		
		# Map contigs to taxid
		taxmap_df = pd.read_csv(tax_map, sep='\t', header=0)
		taxid2name_df = pd.read_csv(taxid_file, sep='\t', header=0)
		taxid2name_df['sp_taxid'] = [int(x) for x in taxid2name_df['@@TAXID']]
		taxid2name_df['sp_name'] = [x.split('|')[-1] for x in taxid2name_df['TAXPATHSN']]

		# Save Predicted SAG contigs to a fasta file
		fasta_out_file = join(save_path, sag_id + '.predicted_contigs.fasta')
		taxa_cnt_dict = {}
		with open(fasta_out_file, 'w') as fasta_out:
			for header, seq in zip(mg_headers, mg_subs):
				contig_id = header.rsplit('|', 1)[0]
				tax_id = taxmap_df.loc[taxmap_df['@@SEQUENCEID'] == contig_id][
																		'TAXID'].values[0]
				if header in SAG_pred_list:
					new_header = '>' + header + '|' + str(tax_id)
					fasta_out.write('\n'.join([new_header, seq]) + '\n')
					if tax_id in taxa_cnt_dict.keys():
						taxa_cnt_dict[tax_id] += 1
					else:
						taxa_cnt_dict[tax_id] = 1
		print('[SAG+]: Predicted subcontigs saved to %s' % basename(fasta_out_file))
		taxa_cnt_df = pd.DataFrame.from_dict(taxa_cnt_dict, orient='index',
												columns=['count']
												)
		taxa_cnt_df['sp_name'] = [list(taxid2name_df.loc[taxid2name_df['sp_taxid'] == int(x)
									]['sp_name'].values)[0] for x in taxa_cnt_df.index
									]
		percent_list = []
		for k in taxa_cnt_dict.keys():
			v = taxa_cnt_dict[k]
			percent_taxa = str(round((v/len(SAG_pred_list))*100, 2))
			percent_list.append(percent_taxa)
			print('[SAG+]: %s makes up %s percent of predicted SAG+' % (k, percent_taxa))
		taxa_cnt_df['percent_recruit'] = percent_list
		taxa_cnt_df.sort_values(by=['count'], inplace=True, ascending=False)
		taxa_tracking_list.append(taxa_cnt_df)
		print('[SAG+]: Completed analysis of %s and %s' % (sag_id, mg_id))

	taxa_tracking_df = pd.concat(taxa_tracking_list)
	taxa_tracking_df.to_csv(join(save_path, 'total_recruit_stats.tsv'), sep='\t')	
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