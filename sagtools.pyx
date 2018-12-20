from itertools import product
import pandas as pd

#  TODO: bring in other defs, imports, cythonize

class seqman:
	def __init__(seq_headers, contig_subseqs, comp_hash_set, kmer_L,
					seq_list, n, o_lap):
		self.seq_headers = seq_headers
		self.contig_subseqs = contig_subseqs
		self.comp_hash_set = comp_hash_set
		self.kmer_L = kmer_L
		self.seq_list = seq_list
		self.l_max = l_max
		self.o_lap = o_lap
		self.seq = seq

		
	def kmer_ID_filter(self.seq_headers, self.contig_subseqs, self.comp_hash_set,
						self.kmer_L
						):
		"Performs kmer hash comparison between comp_hash_set and hashes from contig_subseqs"
		cdef int 
		pass_list = []
		for header, frag in zip(self.seq_headers, self.contig_subseqs):
			tmp, Ls = kmer_slide([(header, frag)], self.kmer_L, self.kmer_L - 1)
			hashes = calc_seg(Ls)
			hashes.sort(reverse=True)
			hashes_set = set(hashes)
			if self.comp_hash_set.intersection(hashes_set):
				pass_list.append(header)

		return pass_list


	def kmer_slide(self.seq_list, self.l_max, self.o_lap):
		"Builds kmers of len l_max and overlap o_lap"
		all_sub_seqs = []
		all_sub_headers = []
		for seq_tup in self.seq_list:
			header, seq = seq_tup
			clean_seq = seq.strip('\n').lower()
			sub_list = get_frags(clean_seq, self.l_max, self.o_lap)
			sub_headers = [header + '_' + str(i) for i, x in enumerate(sub_list, start=0)]
			all_sub_seqs.extend(sub_list)
			all_sub_headers.extend(sub_headers)

		return all_sub_headers, all_sub_seqs


	def get_frags(self.seq, self.l_max, self.o_lap):
		"Fragments seq into subseqs of length l_max and overlap of o_lap"
		"Leftover tail overlaps with tail-1"
		seq_frags = []
		if (self.l_max != 0) and (len(self.seq) > self.l_max):
			offset = self.l_max - self.o_lap
			for i in range(0, len(self.seq), offset):
				if i + self.l_max < len(self.seq):
					frag = self.seq[i:i + self.l_max]
				else:
					frag = self.seq[-self.l_max:]
				seq_frags.append(frag)
		else:
			seq_frags.append(self.seq)

		return seq_frags


	def tetra_cnt(self.seq_list):
		"Build canonical tetranucleotide frequency dataframe from seq list"
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