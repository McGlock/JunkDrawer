import sys
from os import listdir, makedirs, path
from os.path import isfile, join, isdir, basename, dirname
import sourmash
from Bio import SeqIO
import pandas as pd


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
	mg_rpkm_file = '/home/rmclaughlin/Ryan/CAMI_gold/CAMI_I_HIGH/CAMI_high_GoldStandardAssembly.rpkm.tsv'
	max_contig_len = 10000
	overlap_len = 2000
	save_path = '/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/'
	mocksag_path = join(save_path, 'mockSAGs')
	subcontig_path = join(save_path, 'subcontigs')
	sig_path = join(save_path, 'signatures')
	mhr_path = join(save_path, 'minhash_recruits')
	ara_path = join(save_path, 'rpkm_recruits')
	tra_path = join(save_path, 'tetra_recruits')
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


	# Find the SAGs!
	if isdir(sag_path):
		print('[SAG+]: Directory specified, looking for SAGs')
		sag_list = [join(sag_path, f) for f in
					listdir(sag_path) if (f.split('.')[-1] == 'fasta' or f.split('.')[-1] == 'fna')
					]
	elif isfile(sag_path):
		print('[SAG+]: File specified, processing %s' % basename(sag_path))
		sag_list = [sag_path]

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
	else:
		print('[SAG+]: Building subcontigs for %s' % mg_id)
		mg_contigs = get_seqs(mg_file)
		mg_headers, mg_subs = kmer_slide(sag_contigs, max_contig_len,
											overlap_len
											)
		mg_sub_tup = zip(mg_headers, mg_subs)
		with open(join(subcontig_path, mg_id + '.subcontigs.fasta'), 'w') as sub_out:
				sub_rec_list = ['\n'.join(['>'+rec[0], rec[1]]) for rec in mg_sub_tup]
				sub_out.write('\n'.join(sub_rec_list))
	

	#####################################################################################
	###########################                               ###########################
	########################### MinHash Recruitment Algorithm ###########################
	###########################                               ###########################
	#####################################################################################
	# TODO: only load sigs if the recruit files are missing
	# TODO: remove sag_id column in save file
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
	sag_pass_dict = {}
	for sag_id, sag_sub_tup in sag_contigs_dict.items():
		if isfile(join(mhr_path, sag_id + '.mhr_recruits.tsv')):
			print('[SAG+]: Loading  %s and MetaG signature recruit list' % sag_id)
			with open(join(mhr_path, sag_id + '.mhr_recruits.tsv'), 'r') as mhr_in:
				pass_list = [x.split('\t') for x in mhr_in.readlines()]
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
			sag_sig = sag_sig_dict[sag_id]
			for j, mg_sig in enumerate(mg_sig_list):
				jacc_sim = mg_sig.contained_by(sag_sig)
				if jacc_sim >= 1.0:
					pass_list.append((sag_id, mg_sig.name()))
			with open(join(mhr_path, sag_id + '.mhr_recruits.tsv'), 'w') as mhr_out:
				mhr_out.write('\n'.join(['\t'.join(x) for x in pass_list]))
		contig_pass_list = [x[1].rsplit('_', 1)[0] for x in pass_list]
		sag_pass_dict[sag_id] = contig_pass_list
		print('[SAG+]: Recruited %s subcontigs to %s' % (len(pass_list), sag_id))

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
	# TODO: rebuild this algorithm to accept standard output from 'rpkm'
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

	# get MinHash "passed" mg rpkms
	rpkm_stats_list = []
	for sag_id, contig_pass_list in sag_pass_dict.items():
		print('[SAG+]: Calulating/Loading RPKM stats for %s' % sag_id)
		mg_rpkm_pass_df = mg_rpkm_trim_df[mg_rpkm_trim_df.index.isin(contig_pass_list)]
		mg_rpkm_pass_stat_df = mg_rpkm_pass_df.mean().reset_index()
		mg_rpkm_pass_stat_df.columns = ['sample_id', 'mean']
		mg_rpkm_pass_stat_df['std'] = list(mg_rpkm_pass_df.std())
		mg_rpkm_pass_stat_df['var'] = list(mg_rpkm_pass_df.var())
		mg_rpkm_pass_stat_df['IQ_25'] = list(mg_rpkm_pass_df.quantile(0.25))
		mg_rpkm_pass_stat_df['IQ_75'] = list(mg_rpkm_pass_df.quantile(0.75))
		mg_rpkm_pass_stat_df['IQ_10'] = list(mg_rpkm_pass_df.quantile(0.10))
		mg_rpkm_pass_stat_df['IQ_90'] = list(mg_rpkm_pass_df.quantile(0.90))
		# Use passed MG from MHR to recruit more seqs
		iqr_pass_df = mg_rpkm_trim_df.copy()
		for i, col_nm in enumerate(mg_rpkm_trim_df.columns):
			pass_stats = mg_rpkm_pass_stat_df.iloc[[i]]
			pass_min = pass_stats['IQ_10'].values[0]
			pass_max = pass_stats['IQ_90'].values[0]
			iqr_pass_df = iqr_pass_df.loc[(iqr_pass_df[col_nm] >= pass_min) &
											(iqr_pass_df[col_nm] <= pass_max)
											]
		contig_pass_list.extend(iqr_pass_df.index.values)
		iqr_set_list = list(set(contig_pass_list))
		print('[SAG+]: Recruited %s subcontigs to %s' %
				(len(iqr_set_list), sag_id)
				)
		with open(join(ara_path, sag_id + '.ara_recruits.tsv'), 'w') as ara_out:
			ara_out.write('\n'.join(iqr_set_list))

	#####################################################################################
	#####################################################################################
	#####################################################################################
	#####################################################################################


	#####################################################################################
	##################                                                 ##################
	################## Tetranucleotide Frequency Recruitment Algorithm ##################
	##################                                                 ##################
	#####################################################################################

	# Build/Load tetramers for SAGs and MG subset by ara recruits
	for sag_id, sag_sub_tup in sag_subcontigs_dict.items():
		sag_headers = sag_sub_tup[0]
		sag_subs = sag_sub_tup[1]
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

		if isfile(join(tra_path, mg_id + '.' + sag_id + '.tetras.tsv')):
			print('[SAG+]: Loading tetramer Hz matrix for %s' % mg_id + '.' + sag_id)
			mg_tetra_df = pd.read_csv(join(tra_path, mg_id + '.' + sag_id + '.tetras.tsv'),
										sep='\t',index_col=0, header=0
										)
		else:
			print('[SAG+]: Calculating tetramer Hz matrix for %s' % mg_id + '.' + sag_id)
			mg_ara_headers, mg_ara_subs = zip(*[x for x in mg_sub_tup
									if iqr_rpkm_keep_dict[x[0].rsplit('_', 1)[0]] == True
												])
			mg_tetra_df = pd.DataFrame.from_dict(tetra_cnt(mg_ara_subs))
			mg_tetra_df['contig_id'] = mg_ara_headers
			mg_tetra_df.set_index('contig_id', inplace=True)
			mg_tetra_df.to_csv(join(tra_path, mg_id + '.' + sag_id + '.tetras.tsv'),
								sep='\t'
								)



if __name__ == "__main__":
	main()

