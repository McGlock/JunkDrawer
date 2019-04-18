import pandas as pd
from functools import reduce


csv_list = ['RH_S001__insert_270.rpkm.csv', 'RH_S002__insert_270.rpkm.csv',
			'RH_S003__insert_270.rpkm.csv', 'RH_S004__insert_270.rpkm.csv',
			'RH_S005__insert_270.rpkm.csv']

rpkm_df_list = []
for csv in csv_list:
	rpkm_df = pd.read_csv(csv, sep=',', header=0)
	rpkm_df.set_index('Sequence_name', inplace=True)
	sample_id = str(csv.split('_')[1][3])
	rpkm_df.columns = [x + '_' + sample_id for x in rpkm_df.columns]
	rpkm_df_list.append(rpkm_df)

merged_df = reduce(lambda  left,right: pd.merge(left,right,on=['Sequence_name'],
                                            how='outer'), rpkm_df_list)
merged_df.to_csv('CAMI_high_GoldStandardAssembly.rpkm.tsv', sep='\t',
					header=True, index=True
					)

