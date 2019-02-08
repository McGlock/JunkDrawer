import pandas as pd
from functools import reduce


csv_list = ['1021_AC_run134.final.scaffolds.gt1kb.S001.rpkm.csv',
			'1021_AC_run134.final.scaffolds.gt1kb.S002.rpkm.csv',
			'1021_AC_run134.final.scaffolds.gt1kb.S003.rpkm.csv',
			'1021_AC_run134.final.scaffolds.gt1kb.S004.rpkm.csv',
			'1021_AC_run134.final.scaffolds.gt1kb.S005.rpkm.csv']

rpkm_df_list = []
for csv in csv_list:
	rpkm_df = pd.read_csv(csv, sep=',', header=0)
	rpkm_df.set_index('Sequence_name', inplace=True)
	sample_id = str(csv.split('.')[-3][-1])
	rpkm_df.columns = [x + '_' + sample_id for x in rpkm_df.columns]
	rpkm_df_list.append(rpkm_df)

merged_df = reduce(lambda  left,right: pd.merge(left,right,on=['Sequence_name'],
                                            how='outer'), rpkm_df_list)
merged_df.to_csv('1021_AC_run134.final.scaffolds.gt1kb.rpkm.tsv', sep='\t',
					header=True, index=True
					)

