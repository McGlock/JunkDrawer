import sys
from os import listdir, makedirs, path
from os.path import isfile, join, isdir, basename, dirname
import sourmash
from Bio import SeqIO
import pandas as pd
from itertools import product, islice
import umap
from sklearn.mixture import GaussianMixture as GMM
from sklearn.preprocessing import normalize
import numpy as np
from collections import Counter
from subprocess import Popen, PIPE
from sklearn.mixture import BayesianGaussianMixture as BayGMM




def kmer_slide(seq_list, n, o_lap):
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
			seq_frags.append(frag)
	else:
		seq_frags.append(seq)

	return seq_frags


def get_seqs(fasta_file):
	sag_contigs = []
	with open(fasta_file, 'r') as fasta_in:
		for record in SeqIO.parse(fasta_in, 'fasta'):
			f_id = record.id
			f_description = record.description
			f_seq = str(record.seq)
			if f_seq != '':
				sag_contigs.append((f_id, f_seq))

	return sag_contigs


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


def tetra_cnt(seq_list):
	# Dict of all tetramers
	tetra_cnt_dict = {''.join(x):[] for x in product('atgc', repeat=4)}
	# count up all tetramers and also populate the tetra dict
	for seq in seq_list:
		tmp_dict = {k: 0 for k, v in tetra_cnt_dict.items()}
		clean_seq = seq.strip('\n').lower()
		kmer_list = [''.join(x) for x in get_kmer(clean_seq, 4)]
		tetra_counter = Counter(kmer_list)
		total_kmer_cnt = sum(tetra_counter.values())
		# add counter to tmp_dict
		for tetra in tmp_dict.keys():
			count_tetra = int(tetra_counter[tetra])
			tmp_dict[tetra] = count_tetra
		# map tetras to their reverse tetras (not compliment)
		dedup_dict = {}
		for tetra in tmp_dict.keys():
			if (tetra not in dedup_dict.keys()) & (tetra[::-1]
				not in dedup_dict.keys()
				):
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


def mock_SAG(fasta_file):
	# currently just returns half of the genome as a mock SAG
	genome_contigs = get_seqs(fasta_file)
	if len(genome_contigs) != 1:
		half_list = genome_contigs[::2]
	else:
		header = genome_contigs[0][0]
		seq = genome_contigs[0][1]
		half_list = [(header,seq[:int(len(seq)/2)])]

	return half_list


def main():

	sag_path = '/home/rmclaughlin/Ryan/CAMI_gold/CAMI_I_HIGH/source_genomes/'
	mg_file = '/home/rmclaughlin/Ryan/CAMI_gold/CAMI_I_HIGH/CAMI_high_GoldStandardAssembly.fasta'
	mg_rpkm_file = '/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/RPKMs/CAMI_high_GoldStandardAssembly.rpkm.tsv'
	max_contig_len = 10000
	overlap_len = 2000
	save_path = '/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/'
	mocksag_path = join(save_path, 'mockSAGs')
	subcontig_path = join(save_path, 'subcontigs')
	sig_path = join(save_path, 'signatures')
	mhr_path = join(save_path, 'minhash_recruits')
	ara_path = join(save_path, 'rpkm_recruits')
	tra_path = join(save_path, 'tetra_recruits')
	final_path = join(save_path, 'final_recruits')
	ext_path = join(save_path, 'extend_SAGs')
	asm_path = join(save_path, 're-assembled')
	check_path = join(save_path, 'checkM')



	contig_tax_map = '/home/rmclaughlin/Ryan/CAMI_gold/CAMI_I_HIGH/gsa_mapping_pool.binning.trimmed'
	sag_tax_map = '/home/rmclaughlin/Ryan/CAMI_gold/CAMI_I_HIGH/genome_taxa_info.tsv'
	num_components = 20

	# Check if dirs exist, make them if they don't
	if not path.exists(save_path):
		makedirs(save_path)
	if not path.exists(mocksag_path):
		makedirs(mocksag_path)
	if not path.exists(subcontig_path):
		makedirs(subcontig_path)
	if not path.exists(sig_path):
		makedirs(sig_path)
	if not path.exists(mhr_path):
		makedirs(mhr_path)
	if not path.exists(ara_path):
		makedirs(ara_path)
	if not path.exists(tra_path):
		makedirs(tra_path)
	if not path.exists(final_path):
		makedirs(final_path)
	if not path.exists(ext_path):
		makedirs(ext_path)
	if not path.exists(asm_path):
		makedirs(asm_path)
	if not path.exists(check_path):
		makedirs(check_path)


	# Find the SAGs!
	if isdir(sag_path):
		print('[SAG+]: Directory specified, looking for SAGs')
		sag_list = [join(sag_path, f) for f in
					listdir(sag_path) if (f.split('.')[-1] == 'fasta' or
					f.split('.')[-1] == 'fna')
					]
	elif isfile(sag_path):
		print('[SAG+]: File specified, processing %s' % basename(sag_path))
		sag_list = [sag_path]
	# TODO: subcontig function has issue with trailing kmers, needs to stop at last 10K
	# Build Mock SAGs (for testing only), else extract all SAG contigs and headers
	test = True
	print('[SAG+]: Loading/Building subcontigs for all SAGs')
	sag_contigs_dict = {}
	sag_subcontigs_dict = {}
	for sag_file in sag_list:
		sag_basename = basename(sag_file)
		sag_id = sag_basename.rsplit('.', 1)[0]
		if test == True: # (for testing only)
			if isfile(join(mocksag_path, sag_id + '.mockSAG.fasta')):
				sag_contigs = get_seqs(join(mocksag_path, sag_id + '.mockSAG.fasta'))
			else:
				sag_contigs = mock_SAG(sag_file)
				with open(join(mocksag_path, sag_id + '.mockSAG.fasta'), 'w') as mock_out:
					seq_rec_list = ['\n'.join(['>'+rec[0], rec[1]]) for rec in sag_contigs]
					mock_out.write('\n'.join(seq_rec_list))
		else:
			sag_contigs = get_seqs(fasta_file)
		sag_contigs_dict[sag_id] = sag_contigs

		# Build sub sequences for each SAG contig
		if isfile(join(subcontig_path, sag_id + '.subcontigs.fasta')):
			sag_headers, sag_subs = zip(*get_seqs(
										join(subcontig_path, sag_id + '.subcontigs.fasta')
										))
		else:
			sag_headers, sag_subs = kmer_slide(sag_contigs, max_contig_len,
												overlap_len
												)
			with open(join(subcontig_path, sag_id + '.subcontigs.fasta'), 'w') as sub_out:
					sub_rec_list = ['\n'.join(['>'+rec[0], rec[1]]) for rec in
										zip(sag_headers, sag_subs)
										]
					sub_out.write('\n'.join(sub_rec_list))
		sag_subcontigs_dict[sag_id] = sag_headers, sag_subs
	
	# Build/Load subcontigs for Metagenome
	mg_basename = basename(mg_file)
	mg_id = mg_basename.split('.')[0]
	if isfile(join(subcontig_path, mg_id + '.subcontigs.fasta')):
		print('[SAG+]: Loading subcontigs for %s' % mg_id)
		mg_headers, mg_subs = zip(*get_seqs(
									join(subcontig_path, mg_id + '.subcontigs.fasta')
									))
		mg_sub_tup = list(zip(mg_headers, mg_subs))

	else:
		print('[SAG+]: Building subcontigs for %s' % mg_id)
		mg_contigs = get_seqs(mg_file)
		mg_headers, mg_subs = kmer_slide(mg_contigs, max_contig_len,
											overlap_len
											)
		mg_sub_tup = list(zip(mg_headers, mg_subs))
		with open(join(subcontig_path, mg_id + '.subcontigs.fasta'), 'w') as sub_out:
				sub_rec_list = ['\n'.join(['>'+rec[0], rec[1]]) for rec in mg_sub_tup]
				sub_out.write('\n'.join(sub_rec_list))
	

	#####################################################################################
	###########################                               ###########################
	########################### MinHash Recruitment Algorithm ###########################
	###########################                               ###########################
	#####################################################################################
	print('[SAG+]: Starting MinHash Recruitment Algorithm')

	# Calculate/Load MinHash Signatures with SourMash for MG subseqs
	if isfile(join(sig_path, mg_id + '.metaG.sig')): 
		print('[SAG+]: Loading %s Signatures' % mg_id)
		mg_sig_list = sourmash.signature.load_signatures(join(sig_path, mg_id + \
															'.metaG.sig')
															)
	else:
		print('[SAG+]: Building Signatures for %s' % mg_id)
		mg_sig_list = []
		for mg_head, seq in mg_sub_tup:
			up_seq = seq.upper()
			mg_minhash = sourmash.MinHash(n=0, ksize=51, scaled=100)
			mg_minhash.add_sequence(up_seq, force=True)
			mg_sig = sourmash.SourmashSignature(mg_minhash, name=mg_head)
			mg_sig_list.append(mg_sig)
		with open(join(sig_path, mg_id + '.metaG.sig'), 'w') as mg_out:
			sourmash.signature.save_signatures(mg_sig_list,	fp=mg_out)

	# Load comparisons OR Compare SAG sigs to MG sigs to find containment
	print('[SAG+]: Comparing Signatures of SAGs to MetaG contigs')
	minhash_pass_list = []
	for sag_id, sag_sub_tup in sag_contigs_dict.items():
		if isfile(join(mhr_path, sag_id + '.mhr_recruits.tsv')):
			print('[SAG+]: Loading  %s and MetaG signature recruit list' % sag_id)
			with open(join(mhr_path, sag_id + '.mhr_recruits.tsv'), 'r') as mhr_in:
				pass_list = [x.rstrip('\n').split('\t') for x in mhr_in.readlines()]
		else:
			# Calculate\Load MinHash Signatures with SourMash for SAG subseqs
			if isfile(join(sig_path, sag_id + '.SAG.sig')):
				print('[SAG+]: Loading Signature for %s' % sag_id)
				sag_sig = sourmash.signature.load_one_signature(join(sig_path,
																	sag_id + '.SAG.sig')
																	)
			else:
				print('[SAG+]: Building Signature for %s' % sag_id)
				sag_minhash = sourmash.MinHash(n=0, ksize=51, scaled=100)
				for sag_head, sag_subseq in sag_sub_tup:
					sag_upseq = sag_subseq.upper()
					sag_minhash.add_sequence(sag_upseq, force=True)
				sag_sig = sourmash.SourmashSignature(sag_minhash, name=sag_id)
				with open(join(sig_path, sag_id + '.SAG.sig'), 'w') as sags_out:
					sourmash.signature.save_signatures([sag_sig], fp=sags_out)
			print('[SAG+]: Comparing  %s and MetaG signatures' % sag_id)
			pass_list = []
			mg_sig_list = list(mg_sig_list)
			for j, mg_sig in enumerate(mg_sig_list):
				jacc_sim = mg_sig.contained_by(sag_sig)
				mg_nm = mg_sig.name()
				if jacc_sim >= 1.0:
					pass_list.append([sag_id, mg_nm, mg_nm.rsplit('_', 1)[0]])

			with open(join(mhr_path, sag_id + '.mhr_recruits.tsv'), 'w') as mhr_out:
				mhr_out.write('\n'.join(['\t'.join(x) for x in pass_list]))
		minhash_pass_list.extend(pass_list)
		print('[SAG+]: Recruited subcontigs to %s' % sag_id)

	minhash_df = pd.DataFrame(minhash_pass_list, columns=['sag_id', 'subcontig_id',
															'contig_id'
															])

	#####################################################################################
	#####################################################################################
	#####################################################################################
	#####################################################################################


	#####################################################################################
	##########################                                 ##########################
	########################## Abundance Recruitment Algorithm ##########################
	##########################                                 ##########################
	#####################################################################################
	# NOTE: This is built to accept output from 'join_rpkm_out.py' script
	# TODO: Add RPKM cmd call to run within this script
	# TODO: OR impliment Salmon TPM calculation software?
	# TODO: Limit MinHash pass list to contigs with >= 0.51 subcontigs recruited
	print('[SAG+]: Starting Abundance Recruitment Algorithm')

	print('[SAG+]: Loading RPKM values for %s' % mg_id)
	mg_rpkm_df = pd.read_csv(mg_rpkm_file, sep='\t', header=0)
	mg_rpkm_col_list = ['Sequence_name']
	for col in mg_rpkm_df.columns:
		if 'RPKM' in col:
			mg_rpkm_col_list.append(col)
	mg_rpkm_trim_df = mg_rpkm_df[mg_rpkm_col_list]
	mg_rpkm_trim_df = mg_rpkm_trim_df.loc[mg_rpkm_trim_df['Sequence_name']
											!= 'UNMAPPED'
											]
	mg_rpkm_trim_df.set_index('Sequence_name', inplace=True)
	# Normalize data
	normed_rpkm_df = pd.DataFrame(normalize(mg_rpkm_trim_df.values),
								columns=mg_rpkm_trim_df.columns,
								index=mg_rpkm_trim_df.index
								)

	# get MinHash "passed" mg rpkms
	rpkm_pass_list = []
	for sag_id in set(minhash_df['sag_id']):
		print('[SAG+]: Calulating/Loading RPKM stats for %s' % sag_id)
		if isfile(join(ara_path, sag_id + '.ara_recruits.tsv')):
			with open(join(ara_path, sag_id + '.ara_recruits.tsv'), 'r') as ara_in:
				pass_list = [x.rstrip('\n').split('\t') for x in ara_in.readlines()]
		else:
			sag_mh_pass_df = minhash_df[minhash_df['sag_id'] == sag_id]
			mh_cntg_pass_list = set(sag_mh_pass_df['subcontig_id'])
			mg_rpkm_pass_df = normed_rpkm_df[
										normed_rpkm_df.index.isin(mh_cntg_pass_list)
										]
			mg_rpkm_test_df = normed_rpkm_df[
										~normed_rpkm_df.index.isin(mh_cntg_pass_list)
										]
			'''
			n_components = np.arange(1, 100, 1)
			models = [GMM(n, random_state=42)
				  for n in n_components]
			bics = []
			aics = []
			for i, model in enumerate(models):
				n_comp = n_components[i]
				try:
					bic = model.fit(mg_rpkm_pass_df.values,
									mg_rpkm_pass_df.index).bic(mg_rpkm_pass_df.values
									)
					bics.append(bic)
				except:
					print('[WARNING]: BIC failed with %s components' % n_comp)
				try:
					aic = model.fit(mg_rpkm_pass_df.values,
									mg_rpkm_pass_df.index).aic(mg_rpkm_pass_df.values
									)
					aics.append(aic)
				except:
					print('[WARNING]: AIC failed with %s components' % n_comp)

			min_bic_comp = n_components[bics.index(min(bics))]
			min_aic_comp = n_components[aics.index(min(aics))]
			print('[SAG+]: Min AIC/BIC at %s/%s, respectively' % 
					(min_aic_comp, min_bic_comp)
					)

			print('[SAG+]: Using AIC as guide for GMM components')
			print('[SAG+]: Training BayGMM on RPKMs of MinHash recruits')
			bayes_gmm = BayGMM(n_components=min_aic_comp, 
								max_iter=len(mg_rpkm_pass_df.index), random_state=42
								).fit(mg_rpkm_pass_df.values, mg_rpkm_pass_df.index
							)
			sag_predict = bayes_gmm.predict(mg_rpkm_pass_df.values)
			sag_predict_df = pd.DataFrame([x for x in 
											zip(mg_rpkm_pass_df.index, sag_predict)
											], columns=['index', 'label'])
			sag_pred_cnt_df = sag_predict_df.groupby('label').count().reset_index()
			sag_pred_cnt_df.columns = ['label', 'count']
			sag_pred_cnt_df['percent'] = [x/sum(sag_pred_cnt_df['count'])
											for x in sag_pred_cnt_df['count']
											]
			sag_95_labels_df = sag_pred_cnt_df.loc[sag_pred_cnt_df['percent'] > 0.0]
			sag_95_index_list = [x[0] for x in zip(sag_predict_df['index'],
												sag_predict_df['label']
												) if x[1] in list(sag_95_labels_df['label'])
												]
			sag_scores = bayes_gmm.score_samples(mg_rpkm_pass_df.values)
			sag_scores_df = pd.DataFrame(data=sag_scores, index=mg_rpkm_pass_df.index)
			sag_scores_95_df = sag_scores_df.loc[sag_scores_df.index.isin(sag_95_index_list)]
			sag_score_95_min = round(min(sag_scores_95_df.values)[0], 1)
			sag_score_95_max = round(max(sag_scores_95_df.values)[0], 1)
			mg_scores = bayes_gmm.score_samples(mg_rpkm_test_df.values)
			mg_predict = bayes_gmm.predict(mg_rpkm_test_df.values)
			mg_bayes_list = [x for x in zip(mg_scores, mg_predict)]
			mg_bayes_df = pd.DataFrame(data=mg_bayes_list, index=mg_rpkm_test_df.index)
			mg_bayes_df.columns = ['score', 'label']

			if sag_score_95_min == sag_score_95_max:
				print('[SAG+]: Using label instead of score')

				gmm_pass_df = mg_bayes_df.loc[mg_bayes_df['label'
												].isin(sag_95_labels_df['label'])
												]
			else:
				gmm_pass_df = mg_bayes_df.loc[(mg_bayes_df['score'] >= sag_score_95_min) &
												(mg_bayes_df['score'] <= sag_score_95_max) &
												(mg_bayes_df['label'].isin(
																sag_95_labels_df['label']
																))
												]
			
			'''
			'''
			gmm = GMM(n_components=min_aic_comp, random_state=42
							).fit(mg_rpkm_pass_df.values, mg_rpkm_pass_df.index
							)
			sag_scores = gmm.score_samples(mg_rpkm_pass_df.values)
			sag_scores_df = pd.DataFrame(data=sag_scores, index=mg_rpkm_pass_df.index)
			sag_scores_IQR05 = sag_scores_df.quantile(0.05).values[0]
			sag_scores_IQR95 = sag_scores_df.quantile(0.95).values[0]
			sag_score_min = min(sag_scores_df.values)[0]
			sag_score_max = max(sag_scores_df.values)[0]

			mg_scores = gmm.score_samples(mg_rpkm_test_df.values)
			mg_scores_df = pd.DataFrame(data=mg_scores, index=mg_rpkm_test_df.index)
			gmm_pass_df = mg_scores_df.loc[(mg_scores_df[0] >= sag_scores_IQR05) &
											(mg_scores_df[0] <= sag_scores_IQR95)
											]
			'''
			'''
			pass_list = []
			mh_rpkm_indexes = set(list(gmm_pass_df.index.values) + list(mh_cntg_pass_list))
			for md_nm in mh_rpkm_indexes:
				pass_list.append([sag_id, md_nm, md_nm.rsplit('_', 1)[0]])
			print('[SAG+]: Recruited %s subcontigs to %s' % (len(pass_list), sag_id))
			with open(join(ara_path, sag_id + '.ara_recruits.tsv'), 'w') as ara_out:
				ara_out.write('\n'.join(['\t'.join(x) for x in pass_list]))
			<-rpkm_pass_list.extend(pass_list)
			'''
			mg_rpkm_pass_stat_df = mg_rpkm_pass_df.mean().reset_index()
			mg_rpkm_pass_stat_df.columns = ['sample_id', 'mean']
			mg_rpkm_pass_stat_df['std'] = list(mg_rpkm_pass_df.std())
			mg_rpkm_pass_stat_df['var'] = list(mg_rpkm_pass_df.var())
			mg_rpkm_pass_stat_df['skew'] = list(mg_rpkm_pass_df.skew())
			mg_rpkm_pass_stat_df['kurt'] = list(mg_rpkm_pass_df.kurt())
			mg_rpkm_pass_stat_df['IQ_25'] = list(mg_rpkm_pass_df.quantile(0.25))
			mg_rpkm_pass_stat_df['IQ_75'] = list(mg_rpkm_pass_df.quantile(0.75))
			mg_rpkm_pass_stat_df['IQ_10'] = list(mg_rpkm_pass_df.quantile(0.10))
			mg_rpkm_pass_stat_df['IQ_90'] = list(mg_rpkm_pass_df.quantile(0.90))
			mg_rpkm_pass_stat_df['IQ_05'] = list(mg_rpkm_pass_df.quantile(0.05))
			mg_rpkm_pass_stat_df['IQ_95'] = list(mg_rpkm_pass_df.quantile(0.95))
			mg_rpkm_pass_stat_df['IQ_01'] = list(mg_rpkm_pass_df.quantile(0.01))
			mg_rpkm_pass_stat_df['IQ_99'] = list(mg_rpkm_pass_df.quantile(0.99))
			mg_rpkm_pass_stat_df['IQR'] = mg_rpkm_pass_stat_df['IQ_75'] - \
											mg_rpkm_pass_stat_df['IQ_25']
			mg_rpkm_pass_stat_df['upper_bound'] = mg_rpkm_pass_stat_df['IQ_75'] + \
													(1.5 * mg_rpkm_pass_stat_df['IQR'])
			mg_rpkm_pass_stat_df['lower_bound'] = mg_rpkm_pass_stat_df['IQ_75'] - \
													(1.5 * mg_rpkm_pass_stat_df['IQR'])
			# Use passed MG from MHR to recruit more seqs
			iqr_pass_df = mg_rpkm_test_df.copy()
			for i, col_nm in enumerate(mg_rpkm_test_df.columns):
				pass_stats = mg_rpkm_pass_stat_df.iloc[[i]]
				pass_max = pass_stats['upper_bound'].values[0]
				pass_min = pass_stats['lower_bound'].values[0]
				iqr_pass_df = iqr_pass_df.loc[(iqr_pass_df[col_nm] >= pass_min) &
												(iqr_pass_df[col_nm] <= pass_max)
												]
			pass_list = []
			join_rpkm_recruits = set(list(iqr_pass_df.index) + list(mh_cntg_pass_list))
			for md_nm in join_rpkm_recruits:
				#if ((md_nm in iqr_pass_df.index.values) or (md_nm in mh_cntg_pass_list)):
				pass_list.append([sag_id, md_nm, md_nm.rsplit('_', 1)[0]])
			print('[SAG+]: Recruited %s subcontigs to %s' % (len(pass_list), sag_id))
			
			with open(join(ara_path, sag_id + '.ara_recruits.tsv'), 'w') as ara_out:
				ara_out.write('\n'.join(['\t'.join(x) for x in pass_list]))
		rpkm_pass_list.extend(pass_list)

	rpkm_df = pd.DataFrame(rpkm_pass_list, columns=['sag_id', 'subcontig_id',
													'contig_id'
													])
	
	#####################################################################################
	#####################################################################################
	#####################################################################################
	#####################################################################################


	#####################################################################################
	##################                                                 ##################
	################## Tetranucleotide Frequency Recruitment Algorithm ##################
	##################                                                 ##################
	#####################################################################################
	# TODO: Think about using Minimum Description Length (MDL) instead of AIC/BIC
	#		[Normalized Maximum Likelihood or Fish Information Approximation]
	# TODO: Should SAG be included in training data?
	# Build/Load tetramers for SAGs and MG subset by ara recruits
	
	if isfile(join(tra_path, mg_id + '.tetras.tsv')):
		print('[SAG+]: Loading tetramer Hz matrix for %s' % mg_id)
		mg_tetra_df = pd.read_csv(join(tra_path, mg_id + '.tetras.tsv'),
									sep='\t',index_col=0, header=0
									)
	else:
		print('[SAG+]: Calculating tetramer Hz matrix for %s' % mg_id)
		mg_tetra_df = pd.DataFrame.from_dict(tetra_cnt(mg_subs))
		mg_tetra_df['contig_id'] = mg_headers
		mg_tetra_df.set_index('contig_id', inplace=True)
		mg_tetra_df.to_csv(join(tra_path, mg_id + '.tetras.tsv'),
							sep='\t'
							)
	gmm_pass_list = []
	for sag_id, sag_sub_tup in sag_subcontigs_dict.items():
		sag_headers = sag_sub_tup[0]
		sag_subs = sag_sub_tup[1]

		if isfile(join(tra_path, sag_id + '.tra_recruits.tsv')):
			print('[SAG+]: Loading  %s tetramer Hz recruit list' % sag_id)
			with open(join(tra_path, sag_id + '.tra_recruits.tsv'), 'r') as tra_in:
				pass_list = [x.rstrip('\n').split('\t') for x in tra_in.readlines()]
		else:
			if isfile(join(tra_path, sag_id + '.tetras.tsv')):
				print('[SAG+]: Loading tetramer Hz matrix for %s' % sag_id)
				sag_tetra_df = pd.read_csv(join(tra_path, sag_id + '.tetras.tsv'),
											sep='\t', index_col=0, header=0)
			else:
				print('[SAG+]: Calculating tetramer Hz matrix for %s' % sag_id)
				sag_tetra_df = pd.DataFrame.from_dict(tetra_cnt(sag_subs))
				sag_tetra_df['contig_id'] = sag_headers
				sag_tetra_df.set_index('contig_id', inplace=True)
				sag_tetra_df.to_csv(join(tra_path, sag_id + '.tetras.tsv'), sep='\t')

			# Concat SAGs amd MG for GMM
			mg_rpkm_contig_list = list(rpkm_df.loc[rpkm_df['sag_id'] == sag_id
													]['subcontig_id'].values
													)
			mg_rpkm_pass_index = [x for x in mg_tetra_df.index
							if x in mg_rpkm_contig_list
							]
			mg_rpkm_filter_df = mg_tetra_df.loc[mg_tetra_df.index.isin(mg_rpkm_pass_index)]
			concat_tetra_df = pd.concat([sag_tetra_df, mg_rpkm_filter_df])
			normed_tetra_df = pd.DataFrame(normalize(concat_tetra_df.values),
											columns=concat_tetra_df.columns,
											index=concat_tetra_df.index
											)
			sag_normed_tetra_df = normed_tetra_df[
									normed_tetra_df.index.isin(sag_tetra_df.index)
									]
			mg_normed_tetra_df = normed_tetra_df.loc[
									normed_tetra_df.index.isin(mg_rpkm_filter_df.index)
									]

			# UMAP for Dimension reduction of tetras
			sag_features = sag_normed_tetra_df.values
			sag_targets = sag_normed_tetra_df.index.values
			mg_features = mg_normed_tetra_df.values
			mg_targets = mg_normed_tetra_df.index.values
			normed_features = normed_tetra_df.values
			normed_targets = normed_tetra_df.index.values
			
			print('[SAG+]: Dimension reduction of tetras with UMAP')
			umap_trans = umap.UMAP(n_neighbors=2, min_dist=0.0,
							n_components=num_components, metric='manhattan',
							random_state=42
							).fit_transform(normed_features)

			pc_col_names = ['pc' + str(x) for x in range(1, num_components + 1)]
			umap_df = pd.DataFrame(umap_trans, columns=pc_col_names, index=normed_targets)

			sag_umap_df = umap_df.loc[umap_df.index.isin(sag_tetra_df.index)]
			mg_umap_df = umap_df.loc[umap_df.index.isin(mg_tetra_df.index)]
			n_components = np.arange(1, 100, 1)
			models = [GMM(n, random_state=42)
				  for n in n_components]
			bics = []
			aics = []
			for i, model in enumerate(models):
				n_comp = n_components[i]
				try:
					bic = model.fit(sag_umap_df.values,
									sag_umap_df.index).bic(sag_umap_df.values
									)
					bics.append(bic)
				except:
					print('[WARNING]: BIC failed with %s components' % n_comp)
				try:
					aic = model.fit(sag_umap_df.values,
									sag_umap_df.index).aic(sag_umap_df.values
									)
					aics.append(aic)
				except:
					print('[WARNING]: AIC failed with %s components' % n_comp)

			min_bic_comp = n_components[bics.index(min(bics))]
			min_aic_comp = n_components[aics.index(min(aics))]
			print('[SAG+]: Min AIC/BIC at %s/%s, respectively' % 
					(min_aic_comp, min_bic_comp)
					)
			print('[SAG+]: Using AIC as guide for GMM components')
			print('[SAG+]: Training GMM on SAG tetras')
			gmm = GMM(n_components=min_aic_comp, random_state=42
							).fit(sag_umap_df.values, sag_umap_df.index
							)
			print('[SAG+]: GMM Converged: ', gmm.converged_)
			sag_scores = gmm.score_samples(sag_umap_df.values)
			sag_scores_df = pd.DataFrame(data=sag_scores, index=sag_targets)
			sag_score_min = min(sag_scores_df.values)[0]
			sag_score_max = max(sag_scores_df.values)[0]
			mg_scores = gmm.score_samples(mg_umap_df.values)
			mg_scores_df = pd.DataFrame(data=mg_scores, index=mg_targets)
			gmm_pass_df = mg_scores_df.loc[(mg_scores_df[0] >= sag_score_min) &
											(mg_scores_df[0] <= sag_score_max)
											]
			pass_list = []
			for md_nm in gmm_pass_df.index.values:
				pass_list.append([sag_id, md_nm, md_nm.rsplit('_', 1)[0]])
			print('[SAG+]: Recruited %s subcontigs to %s' % (len(pass_list), sag_id))
			with open(join(tra_path, sag_id + '.tra_recruits.tsv'), 'w') as tra_out:
				tra_out.write('\n'.join(['\t'.join(x) for x in pass_list]))
		gmm_pass_list.extend(pass_list)
	gmm_df = pd.DataFrame(gmm_pass_list, columns=['sag_id', 'subcontig_id', 'contig_id'])

	#####################################################################################
	#####################################################################################
	#####################################################################################
	#####################################################################################


	#####################################################################################
	######################                                         ######################
	###################### Collect the recruits and merge with SAG ######################
	######################                                         ######################
	#####################################################################################
	# TODO: Use full contigs instead of subcontigs for co-asm, reduces asm time for Minimus2? CISA?
	# TODO: check for co-asm files before running
	# Merge MinHash and GMM Tetra (passed first by RPKM)
	mh_gmm_merge_df = minhash_df.merge(gmm_df, how='outer',
										on=['sag_id', 'subcontig_id', 'contig_id']
										)
	mh_gmm_merge_df.to_csv(join(final_path, 'final_recruits.tsv'), sep='\t', index=True)
	'''
	# Count # of subcontigs recruited to each SAG
	mh_gmm_cnt_df = mh_gmm_merge_df.groupby(['sag_id', 'contig_id']).count().reset_index()
	mh_gmm_cnt_df.columns = ['sag_id', 'contig_id', 'subcontig_recruits']
	# Build subcontig count for each MG contig
	mg_contig_list = [x.rsplit('_', 1)[0] for x in mg_headers]
	mg_tot_df = pd.DataFrame(zip(mg_contig_list, mg_headers),
									columns=['contig_id', 'subcontig_id'])
	mg_tot_cnt_df = mg_tot_df.groupby(['contig_id']).count().reset_index()
	mg_tot_cnt_df.columns = ['contig_id', 'subcontig_total']
	mg_recruit_df = mh_gmm_cnt_df.merge(mg_tot_cnt_df, how='left', on='contig_id')
	mg_recruit_df['percent_recruited'] = mg_recruit_df['subcontig_recruits'] / \
											mg_recruit_df['subcontig_total']
	mg_recruit_df.sort_values(by='percent_recruited', ascending=False, inplace=True)
	# Only pass contigs that have the magjority of subcontigs recruited (>= 51%)
	mg_recruit_filter_df = mg_recruit_df.loc[mg_recruit_df['percent_recruited'] >= 0.51]
	
	mg_contig_per_max_df = mg_recruit_filter_df.groupby(['contig_id'])[
											'percent_recruited'].max().reset_index()
	mg_contig_per_max_df.columns = ['contig_id', 'percent_max']
	mg_recruit_max_df = mg_recruit_filter_df.merge(mg_contig_per_max_df, how='left',
													on='contig_id')
	# Now pass contigs that have the maximum recruit % of subcontigs
	mg_max_only_df = mg_recruit_max_df.loc[mg_recruit_max_df['percent_recruited'] >=
											mg_recruit_max_df['percent_max']
											]
	'''
	mg_sub_df = pd.DataFrame(mg_sub_tup, columns=['subcontig_id', 'seq'])
	mg_sub_df['contig_id'] = [x.rsplit('_', 1)[0] for x in mg_sub_df['subcontig_id']]
	for sag_id in set(mh_gmm_merge_df['sag_id']):
		sub_merge_df = mh_gmm_merge_df.loc[mh_gmm_merge_df['sag_id'] == sag_id]
		print('[SAG+]: Recruited %s subcontigs from entire analysis for %s' % 
				(sub_merge_df.shape[0], sag_id)
				)
		with open(join(final_path, sag_id + '.final_recruits.fasta'), 'w') as final_out:
			mg_sub_filter_df = mg_sub_df.loc[mg_sub_df['contig_id'
												].isin(sub_merge_df['contig_id'])
												]
			final_mgsubs_list = ['\n'.join(['>'+x[0], x[1]]) for x in
									zip(mg_sub_filter_df['subcontig_id'],
										mg_sub_filter_df['seq']
										)
									]
			final_out.write('\n'.join(final_mgsubs_list))
		# Combine SAG and final recruits
		with open(join(ext_path, sag_id + '.extend_SAG.fasta'), 'w') as cat_file:
			data = []			
			with open(join(mocksag_path, sag_id + '.mockSAG.fasta'), 'r') as sag_in:
				data.extend(sag_in.readlines())
			with open(join(final_path, sag_id + '.final_recruits.fasta'), 'r') as \
				recruits_in:
				data.extend(recruits_in.readlines())
			join_data = '\n'.join(data).replace('\n\n', '\n')
			cat_file.write(join_data)
		
		# Use SPAdes to co-assemble mSAG and recruits
		print('[SAG+]: Re-assembling SAG with final recruits using SPAdes')
		spades_cmd = ['/home/rmclaughlin/bin/SPAdes-3.13.0-Linux/bin/spades.py',
						'--sc', '-k', '21,33,55,77,99,127', '--careful', '--only-assembler',
						'-o', join(asm_path, sag_id), '--trusted-contigs',
						join(mocksag_path, sag_id + '.mockSAG.fasta'),
						'--s1', join(final_path, sag_id + '.final_recruits.fasta')
						]
		run_spades = Popen(spades_cmd, stdout=PIPE)
		print(run_spades.communicate()[0].decode())
		move_cmd = ['mv', join(join(asm_path, sag_id),'scaffolds.fasta'),
						join(asm_path, sag_id + '.asm.fasta')
						]
						
		run_move = Popen(move_cmd, stdout=PIPE)
		print(run_move.communicate()[0].decode())
		clean_cmd = ['rm', '-rf', join(asm_path, sag_id)]
		run_clean = Popen(clean_cmd, stdout=PIPE)
		print(run_clean.communicate()[0].decode())
		
		'''
		# Use minimus2 to merge the SAG and the recruits into one assembly
		toAmos_cmd = ['/home/rmclaughlin/bin/amos-3.1.0/bin/toAmos', '-s',
						join(ext_path, sag_id + '.extend_SAG.fasta'), '-o',
						join(asm_path, sag_id + '.afg')
						]
		run_toAmos = Popen(toAmos_cmd, stdout=PIPE)
		print(run_toAmos.communicate()[0].decode())
		minimus_cmd = ['/home/rmclaughlin/bin/amos-3.1.0/bin/minimus2',
						join(asm_path, sag_id),
						'-D', 'REFCOUNT=0', '-D', 'OVERLAP=200', '-D', 'MINID=95'
						]
		run_minimus = Popen(minimus_cmd, stdout=PIPE)
		print(run_minimus.communicate()[0].decode())
		clean_cmd = ['rm', '-r', join(asm_path, sag_id + '.runAmos.log'),
						join(asm_path, sag_id + '.afg'),
						join(asm_path, sag_id + '.OVL'),
						join(asm_path, sag_id + '.singletons'),
						join(asm_path, sag_id + '.contig'),
						join(asm_path, sag_id + '.ovl'),
						join(asm_path, sag_id + '.coords'),
						join(asm_path, sag_id + '.qry.seq'),
						join(asm_path, sag_id + '.delta'),
						join(asm_path, sag_id + '.bnk'),
						join(asm_path, sag_id + '.ref.seq')
						]
		run_clean = Popen(clean_cmd, stdout=PIPE)
		print(run_clean.communicate()[0].decode())
		'''
	
	# Run CheckM on all new rebuilt/updated SAGs
	print('[SAG+]: Checking all new SAG quality using CheckM')
	checkm_cmd = ['checkm', 'lineage_wf', '--tab_table', '-x',
					'fasta', '--threads', '8', '--pplacer_threads', '8', '-f',
					join(check_path, 'checkM_stdout.tsv'), asm_path, check_path
					]
	run_checkm = Popen(checkm_cmd, stdout=PIPE)
	print(run_checkm.communicate()[0].decode())
	
	
if __name__ == "__main__":
	main()

