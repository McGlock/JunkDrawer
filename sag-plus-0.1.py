import sys
from os import listdir, makedirs, path
from os.path import isfile, join, isdir, basename, dirname
import sourmash
from Bio import SeqIO


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
	print('[SAG+]: Loading/Building subcontigs for all SAGs in: %s' % sag_path)
	sag_contigs_dict = {}
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
	
	# Calculate\Load MinHash Signatures with SourMash for SAG subseqs
	if isfile(join(sig_path, 'SAGs.sig')): 
		print('[SAG+]: Loading Signatures for all SAGs in: %s' % sag_path)
		sag_sig_list = sourmash.signature.load_signatures(join(sig_path, 'SAGs.sig'))
	else:
		sag_sig_list = []
		print('[SAG+]: Building Signatures for all SAGs in: %s' % sag_path)
		for sag_id, sag_sub_tup in sag_contigs_dict.items():
			print('[SAG+]: Building Signatures for %s' % sag_id)
			sag_minhash = sourmash.MinHash(n=0, ksize=51, scaled=100)
			for sag_head, sag_subseq in sag_sub_tup:
				sag_upseq = sag_subseq.upper()
				sag_minhash.add_sequence(sag_upseq, force=True)
			sag_sig = sourmash.SourmashSignature(sag_minhash, name=sag_id)
			sag_sig_list.append(sag_sig)
		with open(join(sig_path, 'SAGs.sig'), 'w') as sags_out:
			sourmash.signature.save_signatures(sag_sig_list, fp=sags_out)
	
	# Calculate/Load MinHash Signatures with SourMash for MG subseqs
	if isfile(join(sig_path, mg_id + '.metaG.sig')): 
		print('[SAG+]: Loading %s Signatures' % mg_id)
		mg_sig_list = list(sourmash.signature.load_signatures(join(sig_path, mg_id + \
															'.metaG.sig')
															))
	else:
		print('[SAG+]: Building Signatures for %s' % mg_id)
		mg_sig_list = []
		for mg_head, seq in mg_sub_tup:
			up_seq = seq.upper()
			mg_minhash = sourmash.MinHash(n=0, ksize=51, scaled=100)
			mg_minhash.add_sequence(up_seq, force=True)
			mg_sig = sourmash.SourmashSignature(mg_minhash, name=mg_head)
			mg_sig_list.append(mg_sig)
		with open(join(save_path, mg_id + '.metaG.sig'), 'w') as mg_out:
			sourmash.signature.save_signatures(mg_sig_list,	fp=mg_out)

	# Load comparisons OR Compare SAG sigs to MG sigs to find containment
	print('[SAG+]: Comparing Signatures of SAGs to MetaG contigs')
	sag_pass_dict = {}
	for i, sag_sig in enumerate(sag_sig_list):
		sag_nm = sag_sig.name()
		if isfile(join(mhr_path, sag_nm + '.mhr_recruits.tsv')): 
			with open(join(mhr_path, sag_sig.name() + '.mhr_recruits.tsv'), 'r') as mhr_in:
				pass_list = [x.split('\t') for x in mhr_in.readlines()]
		else:
			pass_list = []
			for j, s2 in enumerate(mg_sig_list):
				mg_sig = mg_sig_list[j]
				jacc_sim = mg_sig.contained_by(sag_sig)
				if jacc_sim >= 1.0:
					pass_list.append((sag_nm, mg_sig.name()))
			with open(join(mhr_path, sag_nm + '.mhr_recruits.tsv'), 'w') as mhr_out:
				mhr_out.write('\n'.join(['\t'.join(x) for x in pass_list]))
		print('[SAG+]: Recruited %s subcontigs to %s' % (len(pass_list), sag_nm))




if __name__ == "__main__":
	main()

