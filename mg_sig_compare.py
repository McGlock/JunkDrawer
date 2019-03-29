import sourmash
import pandas as pd


mg_sigs_file = 'Ryan/SAG-plus/CAMI_I_HIGH/test_sourmash/CAMI_high_GoldStandardAssembly.minhash.sig'
mg_sig_list = sourmash.signature.load_signatures(mg_sigs_file)
mg_sig_list2 = sourmash.signature.load_signatures(mg_sigs_file)
for i, s1 in enumerate(mg_sig_list):
	for j, s2 in enumerate(mg_sig_list2):
		jaccard_sim = s1.jaccard(s2)
		if jaccard_sim > 0.0:
			print(s1.name(), s2.name(), jaccard_sim) 
