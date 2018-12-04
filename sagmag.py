import sys
from itertools import islice


def tetra_freq(seq_list):
	tetra_cnt_dict = {'aaaa': 0, 'aaat': 0, 'aaag': 0, 'aaac': 0, 'aata': 0,
					'aatt': 0, 'aatg': 0, 'aatc': 0, 'aaga': 0, 'aagt': 0,
					'aagg': 0, 'aagc': 0, 'aaca': 0, 'aact': 0, 'aacg': 0,
					'aacc': 0, 'ataa': 0, 'atat': 0, 'atag': 0, 'atac': 0,
					'atta': 0, 'attt': 0, 'attg': 0, 'attc': 0, 'atga': 0,
					'atgt': 0, 'atgg': 0, 'atgc': 0, 'atca': 0, 'atct': 0, 
					'atcg': 0, 'atcc': 0, 'agaa': 0, 'agat': 0, 'agag': 0,
					'agac': 0, 'agta': 0, 'agtt': 0, 'agtg': 0, 'agtc': 0,
					'agga': 0, 'aggt': 0, 'aggg': 0, 'aggc': 0, 'agca': 0,
					'agct': 0, 'agcg': 0, 'agcc': 0, 'acaa': 0, 'acat': 0,
					'acag': 0, 'acac': 0, 'acta': 0, 'actt': 0, 'actg': 0,
					'actc': 0, 'acga': 0, 'acgt': 0, 'acgg': 0, 'acgc': 0, 
					'acca': 0, 'acct': 0, 'accg': 0, 'accc': 0, 'taaa': 0,
					'taat': 0, 'taag': 0, 'taac': 0, 'tata': 0, 'tatt': 0,
					'tatg': 0, 'tatc': 0, 'taga': 0, 'tagt': 0, 'tagg': 0,
					'tagc': 0, 'taca': 0, 'tact': 0, 'tacg': 0, 'tacc': 0,
					'ttaa': 0, 'ttat': 0, 'ttag': 0, 'ttac': 0, 'ttta': 0,
					'tttt': 0, 'tttg': 0, 'tttc': 0, 'ttga': 0, 'ttgt': 0, 
					'ttgg': 0, 'ttgc': 0, 'ttca': 0, 'ttct': 0, 'ttcg': 0,
					'ttcc': 0, 'tgaa': 0, 'tgat': 0, 'tgag': 0, 'tgac': 0,
					'tgta': 0, 'tgtt': 0, 'tgtg': 0, 'tgtc': 0, 'tgga': 0,
					'tggt': 0, 'tggg': 0, 'tggc': 0, 'tgca': 0, 'tgct': 0,
					'tgcg': 0, 'tgcc': 0, 'tcaa': 0, 'tcat': 0, 'tcag': 0,
					'tcac': 0, 'tcta': 0, 'tctt': 0, 'tctg': 0, 'tctc': 0, 
					'tcga': 0, 'tcgt': 0, 'tcgg': 0, 'tcgc': 0, 'tcca': 0,
					'tcct': 0, 'tccg': 0, 'tccc': 0, 'gaaa': 0, 'gaat': 0,
					'gaag': 0, 'gaac': 0, 'gata': 0, 'gatt': 0, 'gatg': 0,
					'gatc': 0, 'gaga': 0, 'gagt': 0, 'gagg': 0, 'gagc': 0,
					'gaca': 0, 'gact': 0, 'gacg': 0, 'gacc': 0, 'gtaa': 0,
					'gtat': 0, 'gtag': 0, 'gtac': 0, 'gtta': 0, 'gttt': 0, 
					'gttg': 0, 'gttc': 0, 'gtga': 0, 'gtgt': 0, 'gtgg': 0,
					'gtgc': 0, 'gtca': 0, 'gtct': 0, 'gtcg': 0, 'gtcc': 0,
					'ggaa': 0, 'ggat': 0, 'ggag': 0, 'ggac': 0, 'ggta': 0,
					'ggtt': 0, 'ggtg': 0, 'ggtc': 0, 'ggga': 0, 'gggt': 0,
					'gggg': 0, 'gggc': 0, 'ggca': 0, 'ggct': 0, 'ggcg': 0,
					'ggcc': 0, 'gcaa': 0, 'gcat': 0, 'gcag': 0, 'gcac': 0, 
					'gcta': 0, 'gctt': 0, 'gctg': 0, 'gctc': 0, 'gcga': 0,
					'gcgt': 0, 'gcgg': 0, 'gcgc': 0, 'gcca': 0, 'gcct': 0,
					'gccg': 0, 'gccc': 0, 'caaa': 0, 'caat': 0, 'caag': 0,
					'caac': 0, 'cata': 0, 'catt': 0, 'catg': 0, 'catc': 0,
					'caga': 0, 'cagt': 0, 'cagg': 0, 'cagc': 0, 'caca': 0,
					'cact': 0, 'cacg': 0, 'cacc': 0, 'ctaa': 0, 'ctat': 0, 
					'ctag': 0, 'ctac': 0, 'ctta': 0, 'cttt': 0, 'cttg': 0,
					'cttc': 0, 'ctga': 0, 'ctgt': 0, 'ctgg': 0, 'ctgc': 0,
					'ctca': 0, 'ctct': 0, 'ctcg': 0, 'ctcc': 0, 'cgaa': 0,
					'cgat': 0, 'cgag': 0, 'cgac': 0, 'cgta': 0, 'cgtt': 0,
					'cgtg': 0, 'cgtc': 0, 'cgga': 0, 'cggt': 0, 'cggg': 0,
					'cggc': 0, 'cgca': 0, 'cgct': 0, 'cgcg': 0, 'cgcc': 0, 
					'ccaa': 0, 'ccat': 0, 'ccag': 0, 'ccac': 0, 'ccta': 0, 
					'cctt': 0, 'cctg': 0, 'cctc': 0, 'ccga': 0, 'ccgt': 0,
					'ccgg': 0, 'ccgc': 0, 'ccca': 0, 'ccct': 0, 'cccg': 0,
					'cccc': 0
					}
	for seq in seq_list:
		kmer_list = get_tetra(seq)
		for kmer in kmer_list:
			print(kmer)
			### START HERE ###
			tetra_cnt_dict[kmer] += tetra_cnt_dict[kmer]
	return tetra_cnt_dict

def get_tetra(seq, n=4):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result



fasta_file = sys.argv[1]

contig_list = []
with open(fasta_file, 'r') as f:
	data = f.readlines()
	for line in data:
		if '>' not in line:
			contig_list.append(line.strip('\n').lower())

fequencies = tetra_freq(contig_list)
print(fequencies)