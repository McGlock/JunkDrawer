import pandas as pd
import os
import sys


# open mapping file for metaG to metaT
map_file = '/home/rmclaughlin/Lulu/ProcessedData/Illumina_metagenomes/TS_output/metaG_to_metaT_mapping.tsv'
map_df = pd.read_csv(map_file, header=0, sep='\t')
map_dict = {x[0]:x[1] for x in zip(map_df['metaT'], map_df['metaG_asm'])}

TPM_dir = '/home/rmclaughlin/Lulu/ProcessedData/PanMAGs/metaT/TPM/'
TPM_list = os.listdir(TPM_dir)
for TPM_sample in TPM_list:
	print(TPM_sample)
	# open TPM output and build DF from representatives
	tpm_path = TPM_dir + '/' + TPM_sample + '/quant.sf'
	rpkm_dict = {}
	drep_dict = {}
	tpm_df = pd.read_csv(tpm_path, header=0, sep='\t')
	tpm_df['rep_header'] = tpm_df['Name']

	# Build mapping DF from representatives to cluster members
	drep_file = '/home/rmclaughlin/Lulu/ProcessedData/PanMAGs/cluster_reps/cluster_reps.mapped.sam'
	drep_list = []
	with open(drep_file, 'r') as r_in: # Have to deal with stupid header characters
		drep_data = r_in.readlines()
		for dline in drep_data:
			if (('@HD' not in dline) & ('@SQ' not in dline) & ('@PG' not in dline)):
				drep_header = dline.split('\t')[0]
				rep_header = dline.split('\t')[2]
				drep_list.append([drep_header, rep_header])

	drep_df = pd.DataFrame(drep_list, columns=['drep_header', 'rep_header'])
	redup_df = pd.merge(drep_df, tpm_df, on='rep_header', how='left')
	redup_df['rep_header'] = redup_df['drep_header']
	redup_df.drop(columns=['drep_header'], inplace=True)
	redup_df = redup_df[['Name', 'Length', 'EffectiveLength', 'TPM', 'NumReads', 'rep_header']]
	derep_rep_df = pd.concat([tpm_df, redup_df])
	derep_rep_df['cluster'] = [x.split('.', 1)[0] for x in derep_rep_df['rep_header']]
	derep_rep_df['rep_header'] = [x.split('.', 1)[1] for x in derep_rep_df['rep_header']]

	magmap_dir = '/home/rmclaughlin/Lulu/ProcessedData/PanMAGs/rep2query_map/'
	magmap_list = os.listdir(magmap_dir)
	magmap_files = [magmap_dir + x for x in magmap_list if '.mapped.sam' in x]
	mag_list = []
	for magmap_file in magmap_files:
		with open(magmap_file, 'r') as r_in: # Have to deal with stupid header characters
			mag_data = r_in.readlines()
			for mline in mag_data:
				if (('@HD' not in mline) & ('@SQ' not in mline) & ('@PG' not in mline)):
					mag_header = mline.split('\t')[0]
					dmag_header = mline.split('\t')[2]
					mag_list.append([mag_header, dmag_header])

	mag_df = pd.DataFrame(mag_list, columns=['mag_header', 'rep_header'])
	mag_derep_df = pd.merge(mag_df, derep_rep_df, on='rep_header', how='left')
	mag_derep_df['rep_header'] = mag_derep_df['mag_header']
	mag_derep_df.drop(columns=['mag_header'], inplace=True)
	mag_derep_df = mag_derep_df[['Name', 'Length', 'EffectiveLength', 'TPM', 'NumReads',
									'rep_header', 'cluster']
									]

	# Concat all DF
	concat_TPM_df = pd.concat([derep_rep_df, mag_derep_df])
	concat_TPM_df['metaG'] = [x.split('|', 1)[0].split('.')[0] for x in
								concat_TPM_df['rep_header']
								]

	# subset the MAG DF with only the paired metaG and metaT
	sub_mag_df = concat_TPM_df.loc[concat_TPM_df['metaG'] == \
									map_dict[TPM_sample]
									]
	sub_mag_df['seq_header'] = [x.split('|', 1)[1] for x in sub_mag_df['rep_header']]
	sub_mag_df['bin'] = [x.split('|', 1)[0].split('.')[-1] for x in sub_mag_df['rep_header']]
	sub_mag_df['metaT'] = TPM_sample

	sub_mag_df = sub_mag_df[['seq_header', 'cluster', 'bin', 'metaG', 'metaT', 'Length',
								'EffectiveLength', 'NumReads', 'TPM'
								]]
	print(sub_mag_df.shape)
	sv_file = '/home/rmclaughlin/Lulu/ProcessedData/PanMAGs/reduped_mag_TPM/' + \
				map_dict[TPM_sample] + '.' + TPM_sample + \
				'.remapped.tsv'
	sub_mag_df.to_csv(sv_file, sep='\t', index=False)



