import pandas as pd
import sys

# rankedlineage.dmp
ranklin_df = pd.read_csv('~/new_taxdump/rankedlineage.dmp',
                           sep='\t\\|\t', names=['tax_id', 'tax_name', 'species',
                           							'genus', 'family', 'order',
                           							'class', 'phylum', 'kingdom',
                           							'superkingdom'
                           							], engine='python'
                           )
ranklin_df['superkingdom'] = [x.replace('\t|', '') for x in ranklin_df['superkingdom']]
ordered_df = ranklin_df[['tax_id', 'tax_name', 'superkingdom', 'phylum', 'class', 'order', 'family',
							'genus', 'species'
							]]
filled_df = ordered_df.fillna('')
filled_df['lineage'] = [';'.join(x).replace(';;', ';').replace(';;', ';').replace(';;', ';')
								for x in zip(filled_df['superkingdom'], 
												filled_df['phylum'],
												filled_df['class'],
												filled_df['order'],
												filled_df['family'],
												filled_df['genus'],
												filled_df['species']
												)]
# typematerial.dmp
typemat_df = pd.read_csv('~/new_taxdump/typematerial.dmp',
                           sep='\t\\|\t', names=['tax_id', 'tax_name', 'type',
                           							'identifier'
                           							], engine='python'
                           	)
typemat_df['identifier'] = [x.replace('\t|', '') for x in typemat_df['identifier']]

# names.dmp
names_df = pd.read_csv('~/new_taxdump/names.dmp',
                           sep='\t\\|\t', names=['tax_id', 'name_txt', 'unique name',
                           							'name class'
                           							], engine='python'
                           	)
names_df['name class'] = [x.replace('\t|', '') for x in names_df['name class']]

name2taxid = {x[0]:x[1] for x in zip(filled_df['tax_name'], filled_df['tax_id'])}
type2taxid = {x[0]:x[1] for x in zip(typemat_df['identifier'], typemat_df['tax_name'])}
next2taxid = {x[0]:x[1] for x in zip(names_df['name_txt'], names_df['tax_id'])}

# input is MP output functional_and_taxonomic_table.txt in results dir
input_file = sys.argv[1]

tax_func_df = pd.read_csv(input_file, sep='\t', header=0)
tax_id_list = []
for x in tax_func_df['taxonomy']:
	try:
		taxid = int(x.rsplit('(', 1)[1].split(')')[0])
	except:
		if x in name2taxid.keys():
			taxid = name2taxid[x]
		elif x in type2taxid.keys():
			taxid = type2taxid[x]
		else:
			taxid = 'none'
			for i,k in enumerate(next2taxid.keys()):
				if x in k:
					taxid = next2taxid[k]
					break
	tax_id_list.append(taxid)

tax_func_df['tax_id'] = tax_id_list

merge_df = pd.merge(tax_func_df, filled_df, how='left', on='tax_id')

trimmed_df = merge_df[['ORF_ID', 'ORF_length', 'start', 'end', 'Contig_Name',
						'Contig_length', 'strand', 'ec', 'taxonomy', 'product',
						'tax_id', 'tax_name', 'lineage'
						]]


trimmed_df.to_csv(input_file.rsplit('.', 1)[0] + '_taxed.tsv', header=True, sep='\t',
					index=False
					)


