#!/usr/bin/env python
import sys


def subset_fq(sub_ls, fq_in, fq_out):
	with open(sub_ls, 'r') as s:
		sub_dat = s.readlines()
	with open(fq_in, 'r') as i:
		data = i.readlines()
	with open(fq_out, 'w') as o:
		for i, line in enumerate(data):
			if '@' in line[0]:
				seq_id = line.strip('@').split(' ')[0]
				if seq_id in sub_dat:
					rec = [line].extend(data[i+1:i+3])
					o.write(rec)


subset_list = sys.argv[1]
input_fq = sys.argv[2]
output_fq = sys.argv[3]
subset_fq(subset_list, input_fq, output_fq)


