import pandas as pd

mag_cluster_file = '../sourmash/MAG_UMAP.tsv'

mag_cluster_df = pd.read_csv(mag_cluster_file, sep='\t', header=0)

cluster_set = set(mag_cluster_df['sm_clusters'])
for cluster in cluster_set:
	print(cluster)
	cluster_df = mag_cluster_df.loc[mag_cluster_df['sm_clusters'] == cluster]
	merge_list = []
	for sample_id in cluster_df['Sample_Bin_ID']:
		ffn_file = './MAG_ffn/' + sample_id + '.ffn'
		with open(ffn_file, 'r') as ffn:
			data = ffn.readlines()
		new_data = []
		for line in data:
			if '>' in line:
				new_line = '>' + sample_id + '|' + \
							line.split(' ', 1)[0].strip('>') + \
							'|' + line.split(' ', 1)[1].replace(' ', '_')
				new_data.append(new_line)
			else:
				new_data.append(line)
		merge_list.extend(new_data)
	cluster_out = './merged_ffn/' + str(cluster) + '.ffn.fasta'
	with open(cluster_out, 'w') as c_out:
		c_out.write(''.join(merge_list))
