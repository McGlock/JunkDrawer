import pandas as pd
import os

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

	rpkm_file = '/home/rmclaughlin/Ryan/Lulu/MAG_RPKMs/metaG/2AD43II_FD.metaG.rpkm.csv'
	rpkm_dict = {}
	drep_dict = {}
	with open(rpkm_file, 'r') as r_in: # Have to deal with stupid header characters
		data = r_in.readlines()
		for line in data[1:]:
			line = line.strip('\n')
			if 'UNMAPPED' not in line:
				query_sample_id = line.split(',', 1)[0].split('/')[-1].split('.', 1)[0]
				seq_type = line.split(',', 1)[0].split('/')[-1].split('.')[1]
				rep_header = line.split(',', 1)[1].rsplit(',', 2)[0]
				rep_cluster_id = line.split(',', 1)[1].rsplit(',', 2)[0].split('|', 1)[0].split('.')[0]
				rep_sample_id = line.split(',', 1)[1].rsplit(',', 2)[0].split('|', 1)[0].split('.')[1]
				rep_sample_id = line.split(',', 1)[1].rsplit(',', 2)[0].split('|', 1)[0].split('.')[1]
				rep_bin_id = line.split(',', 1)[1].rsplit(',', 2)[0].split('|', 1)[0].split('.')[2]
				rep_seq_header = line.split(',', 1)[1].rsplit(',', 2)[0].split('|', 1)[1]
				read_cnt = line.split(',', 1)[1].rsplit(',', 2)[1]
				rpkm_val = line.split(',', 1)[1].rsplit(',', 2)[2]
				rpkm_dict[rep_header] = [query_sample_id, seq_type, rep_header, rep_cluster_id,
											rep_sample_id, rep_bin_id, rep_seq_header, read_cnt,
											rpkm_val
											]
				drep_dict[rep_header] = [query_sample_id, seq_type, rep_header, rep_cluster_id,
											rep_sample_id, rep_bin_id, rep_seq_header, read_cnt,
											rpkm_val
											]

	rep_rpkm_df = pd.DataFrame.from_dict(rpkm_dict, orient='index',
											columns=['query_sample_id', 'seq_type',
													'rep_header', 'rep_cluster_id',
													'rep_sample_id', 'rep_bin_id',
													'rep_seq_header', 'read_cnt',
													'rpkm_val'
													]).reset_index()
	rep_rpkm_df.columns = ['drep_header' if x == 'index' else x for x in rep_rpkm_df.columns]

	drep_file = '/home/rmclaughlin/Ryan/Lulu/MAG_RPKMs/cluster_reps/cluster_reps.mapped.sam'
	drep_list = []
	with open(drep_file, 'r') as r_in: # Have to deal with stupid header characters
		drep_data = r_in.readlines()
		for dline in drep_data:
			if (('@HD' not in dline) & ('@SQ' not in dline) & ('@PG' not in dline)):
				drep_header = dline.split('\t')[0]
				rep_header = dline.split('\t')[2]
				drep_list.append([drep_header] + rpkm_dict[rep_header])
				drep_dict[rep_header] = rpkm_dict[rep_header]
	# Build df from parsed RPKM file
	drep_rpkm_df = pd.DataFrame(drep_list, columns=['drep_header', 'query_sample_id', 'seq_type',
													'rep_header', 'rep_cluster_id',
													'rep_sample_id', 'rep_bin_id',
													'rep_seq_header', 'read_cnt',
													'rpkm_val'
													])


	magmap_dir = '/home/rmclaughlin/Ryan/Lulu/MAG_RPKMs/rep2query_map/'
	magmap_list = os.listdir(magmap_dir)
	magmap_files = [magmap_dir + x for x in magmap_list if '.mapped.sam' in x]
	mag_list = []
	for magmap_file in magmap_files:
		print(magmap_file)
		with open(magmap_file, 'r') as r_in: # Have to deal with stupid header characters
			mag_data = r_in.readlines()
			for mline in mag_data:
				if (('@HD' not in dline) & ('@SQ' not in dline) & ('@PG' not in dline)):
					mag_header = dline.split('\t')[0]
					drep_header = dline.split('\t')[2]
					mag_list.append([mag_header] + drep_dict[drep_header])

	mag_rpkm_df = pd.DataFrame(mag_list, columns=['drep_header', 'query_sample_id', 'seq_type',
													'rep_header', 'rep_cluster_id',
													'rep_sample_id', 'rep_bin_id',
													'rep_seq_header', 'read_cnt',
													'rpkm_val'
													])

	# Now contains all deduped reps and reps that were removed by dedup
	concat_rpkm_df = pd.concat([rep_rpkm_df, drep_rpkm_df, mag_rpkm_df])

	# Trim to only values that are from the same sample from which the bin was built
	trim_rpkm_df = concat_rpkm_df.loc[concat_rpkm_df['query_sample_id'] == 
										concat_rpkm_df['rep_sample_id']
										]
	print(trim_rpkm_df.head())
	print(trim_rpkm_df.shape)
	sv_file = '/home/rmclaughlin/Ryan/Lulu/MAG_RPKMs/reduped_mag_rpkm/' + /
				rep_sample_id + '.' +  seq_type + 'remapped.tsv'
	trim_rpkm_df.to_csv(sv_file, sep='\t', index=False)



