from itertools import product, islice
import pandas as pd

#  TODO: might be static typing to much, should ask someone about that.

class SeqMan:
	"""For messing with sequences"""

	def __init__(self, list seq_headers=None, list contig_subseqs=None,
					set comp_hash_set=None,	int kmer_L=0, list seq_list=None,
					int l_max=0, int o_lap=0, list seq_tup_list=None,
					int n=0, str seq=None
					):
		self.seq_headers = seq_headers
		self.contig_subseqs = contig_subseqs
		self.comp_hash_set = comp_hash_set
		self.kmer_L = kmer_L
		self.seq_list = seq_list
		self.l_max = l_max
		self.o_lap = o_lap
		self.seq = seq
		self.seq_tup_list = seq_tup_list
		self.n = n

		
	def kmer_ID_filter(self, seq_headers, contig_subseqs, comp_hash_set,
						kmer_L
						):
		"""
		Performs kmer hash comparison between comp_hash_set and hashes
		from contig_subseqs
		"""
		cdef list pass_list = []
		cdef str header
		cdef str frag
		cdef list tmp
		cdef list Ls
		cdef list hashes
		cdef set hashes_set

		for header, frag in zip(seq_headers, contig_subseqs):
			tmp, Ls = self.kmer_slide([(header, frag)], kmer_L, kmer_L - 1)
			hashes = self.calc_seg(Ls)
			hashes.sort(reverse=True)
			hashes_set = set(hashes)
			if comp_hash_set.intersection(hashes_set):
				pass_list.append(header)

		return pass_list


	def calc_seg(self, contig_subseqs):
		"""
		Iterates through nucs of a subseq list and sums from calc_nuc
		Creates Perfect hash list for each subseq
		"""
		cdef list seg_list = []
		cdef long long seg_sum
		cdef int i
		cdef str nuc
		cdef list nuc_naughty_list = ['r', 'y', 's', 'w', 'k',
										'm', 'b', 'd', 'h', 'v', 'n']

		for seq in contig_subseqs:
			seg_sum = 0
			for i, nuc in enumerate(seq, start=0):
				if nuc not in nuc_naughty_list:  # TODO: don't know what to do with these
					seg_sum += self.calc_nuc(nuc, i)
			seg_list.append(seg_sum)

		return seg_list


	def calc_nuc(self, nuc, ind):
		"""
		Calculates number for perfect hash from the nucleotide and the location index
		"""
		cdef dict nuc_dict = {'a': 0, 't': 1, 'c': 2, 'g': 3}
		cdef long long nuc_hash

		nuc_hash = nuc_dict[nuc] * (4**ind)

		return nuc_hash


	def kmer_slide(self, seq_tup_list, l_max, o_lap):
		"Builds kmers of len l_max and overlap o_lap"
		cdef list all_sub_seqs = []
		cdef list all_sub_headers = []
		cdef tuple seq_tup
		cdef str header
		cdef str seq
		cdef str clean_seq
		cdef list sub_list
		cdef list sub_headers
		cdef int i
		cdef str x

		for seq_tup in seq_tup_list:
			header, seq = seq_tup
			clean_seq = seq.strip('\n').lower()
			sub_list = self.get_frags(clean_seq, l_max, o_lap)
			sub_headers = [header + '_' + str(i) for i, x in enumerate(sub_list, start=0)]
			all_sub_seqs.extend(sub_list)
			all_sub_headers.extend(sub_headers)

		return all_sub_headers, all_sub_seqs


	def get_frags(self, seq, l_max, o_lap):
		"Fragments seq into subseqs of length l_max and overlap of o_lap"
		"Leftover tail overlaps with tail-1"
		cdef list seq_frags = []
		cdef int offset
		cdef int i
		cdef str frag

		if (l_max != 0) and (len(seq) > l_max):
			offset = l_max - o_lap
			for i in range(0, len(seq), offset):
				if i + l_max < len(seq):
					frag = seq[i:i + l_max]
				else:
					frag = seq[-l_max:]
				seq_frags.append(frag)
		else:
			seq_frags.append(seq)

		return seq_frags


	def tetra_cnt(self, seq_list):
		"Build canonical tetranucleotide frequency dataframe from seq list"
		cdef dict tetra_cnt_dict
		cdef tuple x
		cdef str seq
		cdef dict tmp_dict
		cdef str k
		cdef list v
		cdef int total_kmer_cnt
		cdef str clean_seq
		cdef list kmer_list
		cdef str tetra
		cdef int count_tetra
		cdef dict dedup_dict
		cdef dict tetra_prop_dict
		cdef double t_prop
 

		# Dict of all tetramers
		tetra_cnt_dict = {''.join(x):[] for x in product('atgc', repeat=4)}
		# count up all tetramers and also populate the tetra dict
		for seq in seq_list:
			tmp_dict = {k: 0 for k, v in tetra_cnt_dict.items()}
			total_kmer_cnt = 0
			clean_seq = seq.strip('\n').lower()
			kmer_list = [''.join(x) for x in self.get_kmer(clean_seq, 4)]
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


	def get_kmer(self, seq, n):
		"Returns a sliding window (of width n) over data from the iterable"
		"   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...				   "
		cdef tuple result

		it = iter(seq)
		result = tuple(islice(it, n))
		if len(result) == n:
			yield result
		for elem in it:
			result = result[1:] + (elem,)
			yield result