import sys
from os import listdir, makedirs, path
from os.path import isfile, join, isdir, basename, dirname
import sourmash


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

	sag_path = '~/Ryan/CAMI_gold/CAMI_I_HIGH/source_genomes/'
	mg_file = '~/Ryan/CAMI_gold/CAMI_I_HIGH/CAMI_high_GoldStandardAssembly.fasta'
	mg_rpkm_file = '~/Ryan/CAMI_gold/CAMI_I_HIGH/CAMI_high_GoldStandardAssembly.rpkm.tsv'
	max_contig_len = 10000
	overlap_len = 2000
	save_path = '~/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/'
	contig_tax_map = '~/Ryan/CAMI_gold/CAMI_I_HIGH/gsa_mapping_pool.binning.trimmed'
	sag_tax_map = '~/Ryan/CAMI_gold/CAMI_I_HIGH/genome_taxa_info.tsv'
	num_components = 20

	# Find the SAGs!
	if not path.exists(save_path):
		makedirs(save_path)

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
	sag_contigs_dict = {}
	for sag_file in sag_list:
		sag_basename = basename(sag_file)
		sag_id = sag_basename.rsplit('.', 1)[0]
		if test == True: # (for testing only)
			if isfile(join(save_path, sag_id + '.mockSAG.fasta')): 
				sag_contigs = get_seqs(fasta_file)
			else:
				sag_contigs, sag_raw_contig_headers = mock_SAG(sag_file)
				with open(join(save_path, sag_id + '.mockSAG.fasta'), 'w') as mock_out:
					mock_out
		else:
			sag_contigs = get_seqs(fasta_file)
		sag_contigs_dict[sag_id] = sag_contigs

		# Build sub sequences for each SAG contig
		sag_headers, sag_subs = kmer_slide(sag_contigs, max_contig_len,
												overlap_len
												)
	
	# Build/Load subcontigs for Metagenome
	mg_basename = basename(mg_file)
	mg_id = mg_basename.split('.')[0]
	mg_contigs = get_seqs(mg_file)
	mg_headers, mg_subs = kmer_slide(mg_contigs, max_contig_len,
										overlap_len
										)
	mg_sub_tup = zip(mg_headers, mg_subs)

	# Calculate\Load MinHash Signatures with SourMash for SAG subseqs
	if isfile(join(save_path, 'SAGs.sig')): 
		print('[SAG+]: Loading Signatures for all SAGs in: %s' % sag_path)
		sag_sig_list = sourmash.signature.load_signatures(join(save_path, 'SAGs.sig'))
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
		with open(join(save_path, 'SAGs.sig'), 'w') as sags_out:
			sourmash.signature.save_signatures(sag_sig_list, fp=sags_out)
	
	# Calculate/Load MinHash Signatures with SourMash for MG subseqs
	if isfile(join(save_path, mg_id + '.metaG.sig')): 
		print('[SAG+]: Loading %s Signatures' % mg_id)
		mg_sig_list = list(sourmash.signature.load_signatures(join(save_path, mg_id + \
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

	# Compare SAG sigs to MG sigs to find containment
	print('[SAG+]: Comparing Signatures of SAGs to MetaG contigs')
	sag_pass_dict = {}
	for i, sag_sig in enumerate(sag_sig_list):
		pass_list = []
		for j, s2 in enumerate(mg_sig_list):
			mg_sig = mg_sig_list[j]
			jacc_sim = mg_sig.contained_by(sag_sig)
			if jacc_sim >= 1.0:
				pass_list.append((sag_sig.name(), mg_sig.name()))
		print('[SAG+]: Recruited %s subcontigs to %s' % (len(pass_list), sag_sig.name()))




if __name__ == "__main__":
	main()

