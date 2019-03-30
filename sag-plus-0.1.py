


def main():
	sag_path = sys.argv[1]
	mg_file = sys.argv[2]
	mg_rpkm_file = sys.argv[3]
	max_contig_len = int(sys.argv[4])
	overlap_len = int(sys.argv[5])
	save_path = sys.argv[6]
	contig_tax_map = sys.argv[7]
	sag_tax_map = sys.argv[8]
	num_components = 20 # int(sys.argv[9])

	sq = sagtools.SeqMan()

	# magic numbers
	#num_components = 2

	# TODO: build argv interface
	# TODO: extract all magic numbers to be built into argv
	# TODO: Refactor all definitions, they are a mess :/
	# TODO: output number of tetras/kmers skipped b/c of ambiguous nucs
	# TODO: cythonize tetra_cnt and calc_seg to optimize

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

	##################################################################################
	#################### Refactoring Recruiting Algorithms ###########################
	##################################################################################
	# Build Mock SAGs (for testing only), else extract all SAG contigs and headers
	test = True
	sag_contigs_dict = {}
	for sag_file in sag_list:
		sag_basename = basename(sag_file)
		sag_id = sag_basename.rsplit('.', 1)[0]
		if test == True: # (for testing only)
			sag_contigs, sag_raw_contig_headers = mock_SAG(sag_file)
		else:
			sag_contigs = get_seqs(fasta_file)
		sag_contigs_dict[sag_id] = sag_contigs

		# Build sub sequences for each SAG contig
		sag_headers, sag_subs = sq.kmer_slide(seq_tup_list=sag_contigs,
												l_max=max_contig_len,
												o_lap=overlap_len
												)
	
	# Build/Load subcontigs for Metagenome
	mg_basename = basename(mg_file)
	mg_id = mg_basename.split('.')[0]
	mg_contigs = get_seqs(mg_file)
	mg_headers, mg_subs = sq.kmer_slide(mg_contigs,	max_contig_len,
										overlap_len
										)
	mg_sub_tup = zip(mg_headers, mg_subs)

	# Calculate\Load MinHash Signatures with SourMash for SAG subseqs
	if isfile(join(save_path, 'SAGs.sig')): 
		sag_sig_list = list(sourmash.signature.load_signatures(join(save_path, 'SAGs.sig'))
							)
		print('[SAG+]: Loaded Signatures for all SAGs in: %s' % sag_path)
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
		mg_sig_list = list(sourmash.signature.load_signatures(join(save_path, mg_id + \
															'.metaG.sig')
															))
		print('[SAG+]: Loaded %s Signatures' % mg_id)
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
	sag_pass_dict = {}
	for i, s1 in enumerate(sag_sig_list):
		sag_sig = sag_sig_list[i]
		pass_list = []
		for j, s2 in enumerate(mg_sig_list):
			mg_sig = mg_sig_list[j]
			jacc_sim = mg_sig.contained_by(sag_sig)
			if jacc_sim >= 1.0:
				pass_list.append((sag_sig.name(), mg_sig.name()))
		print('[SAG+]: Identified %s subcontigs with Sourmash' % len(pass_list))
	##################################################################################
	##################################################################################
	##################################################################################
	sys.exit()


if __name__ == "__main__":
	main()

