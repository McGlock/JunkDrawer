
#  TODO: bring in other defs, imports, cythonize

class kmerClass:
    def __init__(seq_headers, contig_subseqs, comp_hash_set, kmer_L):
        self.seq_headers = seq_headers
        self.contig_subseqs = contig_subseqs
        self.comp_hash_set = comp_hash_set
        self.kmer_L = kmer_L

    def kmer_ID_filter(self):
		pass_list = []
		for header, frag in zip(self.seq_headers, self.contig_subseqs):
			tmp, Ls = get_subseqs([(header, frag)], self.kmer_L, self.kmer_L - 1)
			hashes = calc_seg(Ls)
			hashes.sort(reverse=True)
			hashes_set = set(hashes)
			if self.comp_hash_set.intersection(hashes_set):
				pass_list.append(header)

		return pass_list