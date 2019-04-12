import pandas as pd
from functools import reduce


tsv_list = ['RH_S001__insert_270/quant.sf', 'RH_S003__insert_270/quant.sf',
			'RH_S005__insert_270/quant.sf', 'RH_S002__insert_270/quant.sf',
			'RH_S004__insert_270/quant.sf']

tpm_df_list = []
for tsv in tsv_list:
	tpm_df = pd.read_csv(tsv, sep='\t', header=0)
	tpm_df.set_index('Name', inplace=True)
	sample_id = str(tsv.split('/')[0].split('_')[1][-1])
	tpm_df.columns = [x + '_' + sample_id for x in tpm_df.columns]
	tpm_df_list.append(tpm_df)

merged_df = reduce(lambda  left,right: pd.merge(left,right,on=['Name'],
                                            how='outer'), tpm_df_list)
merged_df.to_csv('CAMI_high_GoldStandardAssembly.tpm.tsv', sep='\t',
					header=True, index=True
					)

