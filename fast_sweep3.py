#!/usr/bin/env python
import sys
import pandas as pd


def clean_fasta(fasta_in, fasta_out, tax_id_file):

	tax_id_df = pd.read_csv(tax_id_file, sep='\t', names=['ID', 'spp | acc', 'lineage'])
	tax_id_df['acc'] = [x.split('|')[1] for x in tax_id_df['spp | acc']]
	with open(fasta_in, 'r') as i:
		data = i.readlines()
	with open(fasta_out, 'w') as o:
		for line in data:
			if '>' in line:
				acc = line.strip('>').rstrip('\n')
				ID = tax_id_df[tax_id_df['acc'] == acc].iloc[0]['ID']
				new_line = '>' + str(ID) + '_' + fasta_in.split('/')[-1].split('.')[0] + '\n'
				o.write(new_line)
			else:
				o.write(line)


recruit_fasta_file = sys.argv[1]
clean_fasta_file = sys.argv[2]
tax_id_file = sys.argv[3]
clean_fasta(recruit_fasta_file, clean_fasta_file, tax_id_file)
