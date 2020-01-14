import sys
import pandas as pd



cm_file = sys.argv[1]
cm_df = pd.read_csv(cm_file, sep='\t', header=0)

cm_df['seq_type'] = 'metagenome'
cm_df['header'] = [x.split(' ', 1)[0] for x in cm_df['Query']]

rpkm_file = sys.argv[2]
rpkm_df = pd.read_csv(rpkm_file, sep=',', header=0)
rpkm_df['header'] = rpkm_df['Sequence_name']
rpkm_dict = {x[0]:x[1] for x in zip(rpkm_df['header'], rpkm_df['RPKM'])}
cm_mt_df = cm_df.copy()
cm_mt_df['seq_type'] = 'metatranscriptome'
cm_mt_df['Abundance'] = [rpkm_dict[x] for x in cm_mt_df['header']]

cat_df = pd.concat([cm_df, cm_mt_df])
outfile = cm_file.rsplit('.', 1)[0] + '.metaG_metaT.tsv'
cat_df.to_csv(outfile, sep='\t', index=False)