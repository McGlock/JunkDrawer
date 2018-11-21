#!/usr/bin/env python
import sys


def clean_fasta(fasta_in, fasta_out):
	with open(fasta_in, 'r') as i:
		data = i.readlines()
	with open(fasta_out, 'w') as o:
		for line in data:
			if '>' in line:
				split_line = line.split('/')
				acc = split_line[0]
				if '|' in acc:
					acc = '>' + acc.split('|')[1]
				new_line = acc + ' ' + '_'.join(split_line[1:])
				o.write(new_line)
			else:
				o.write(line)


recruit_fasta_file = sys.argv[1]
clean_fasta_file = sys.argv[2]
clean_fasta(recruit_fasta_file, clean_fasta_file)


