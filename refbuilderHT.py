#!/usr/bin/env python

import sys
from os import listdir, makedirs, chdir, chmod, walk, stat, environ
from os.path import exists, dirname, realpath, join, isfile, isdir
from shutil import copy
from subprocess import Popen
from multiprocessing import Pool 

def makedir(path, dir_name):
	direct = join(path,dir_name)
	if not exists(direct):
		makedirs(direct)
	return direct

def run_cmds(cmds):
	my_env = environ.copy()
	wpath = cmds[0]
	print(wpath)
	err_file = join(wpath, 'stderr.txt')
	out_file = join(wpath, 'stdout.txt')

	with open(err_file,'w') as err:
		with open(out_file, 'w') as out:
			for cmd in cmds[1:]:
				cmd_str = ' '.join(cmd)
				print(cmd_str)
				proc = Popen(cmd_str, stderr=err, stdout=out,
								shell=True, env=my_env)
				proc.communicate()

def hmm_get_len(hmm):
	with open(hmm, 'r') as h:
		data = h.readlines()
		for line in data:
			if 'LENG' in line:
				hmm_l = int(line.split()[-1])
	return hmm_l


# Magic numbers
len_thresh = 0.8 # no seqs < 80% len of the ref hmm

# number of cores on system
num_jobs = sys.argv[4] # multiprocessing.cpu_count()
# set paths
home_path = sys.argv[1]
file_path = sys.argv[2]
refpkg_workdir = sys.argv[3]
ref_db_path = sys.argv[5]
acc2taxid_lineage_file = sys.argv[6]

# create refpkg directory
refpkg_workpath = makedir(home_path, refpkg_workdir)
if isdir(file_path) == True:
	file_list = [join(file_path, f) for f in listdir(file_path) 
					if f.find('.fasta') != -1
					]
	if not file_list:
		file_list = [join(file_path, f) for f in listdir(file_path) 
						if f.find('.hmm') != -1
						]
	if not file_list:
		sys.exit('Warning: Directory is empty...')

elif isfile(file_path) == True:
	if '.fasta' in file_path:
		file_list = [file_path]
	elif '.hmm' in file_path:
		file_list = [file_path]
	else:
		sys.exit('Unsuported filetype for reference data...')



cmds_list = []
# build command list for build refpkgs from ref fasta list
for file in file_list:
	run_name = file.split('/')[-1].split('.')[0]
	run_path = makedir(refpkg_workpath, run_name)
	refpkg_path = makedir(run_path, run_name + '_refpkg')
	ref_aln_file = join(run_path, run_name + '.ref.aln')
	if '.fasta' in file:
		ref_hmm_file = join(run_path, run_name + '.ref.hmm')
	elif '.hmm' in file:
		ref_hmm_file = file
	hmm_len = hmm_get_len(ref_hmm_file)
	len_cutoff = int(hmm_len*len_thresh)
	recruit_aln_file = join(run_path, run_name + '.recruit.sto')
	recruit_fasta_file = join(run_path, run_name + '.recruit.fasta')
	clean_fasta_file = join(run_path, run_name + '.clean.fasta')
	nolin_acc_file = join(run_path, run_name + '.nolin.txt')
	taxed_fasta_file = join(run_path, run_name + '.taxed.fasta')
	sorted_fasta_file = join(run_path, run_name + '.sorted.fasta')
	filtered_fasta_file = join(run_path, run_name + '.filtered.fasta')
	clean2_fasta_file = join(run_path, run_name + '.clean2.fasta')
	cluster_fasta_file = join(run_path, run_name + '.clean2.cluster.fasta')
	tax_id_file = join(run_path, run_name + '.tax_ids.txt')
	clean3_fasta_file = join(run_path, run_name + '.clean3.fasta')
	final_hmm_file = join(refpkg_path, run_name + '.final.hmm')
	final_tree_file = join(refpkg_path, run_name + '.final.tre')
	final_aln_file = join(refpkg_path, run_name + '.final.fa')
	final_tax_id_file = join(refpkg_path, run_name + '.tax_ids.txt')

	# TODO: add all file objects from build_ref_lineage

	align_cmd = ['muscle', '-in', file, '-out', ref_aln_file]
	hmmbuild_cmd = ['hmmbuild --cpu 4', ref_hmm_file, ref_aln_file]
	hmmsearch_cmd = ['hmmsearch --cpu 4 -E 1e-7 -A', recruit_aln_file, ref_hmm_file, ref_db_path]
	sto2fasta_cmd = ['seqmagick convert --ungap', recruit_aln_file,
						recruit_fasta_file
						]
	clean_fasta_cmd = ['python ~/bin/JunkDrawer/fast_sweep.py', recruit_fasta_file,
						clean_fasta_file
						]
	build_ref_lineage_cmd = ['python ~/bin/JunkDrawer/build_ref_lineage.py',
								acc2taxid_lineage_file, clean_fasta_file, run_path
								]
	sort_fasta_cmd = ['seqmagick convert --sort length-desc', taxed_fasta_file,
						sorted_fasta_file
						]
	filter_fasta_cmd = ['seqmagick convert', sorted_fasta_file, filtered_fasta_file,
						'--exclude-from-file', nolin_acc_file, '--min-ungapped-length',
						str(len_cutoff), '--deduplicate-taxa', '--deduplicate-sequences'
						]
	clean2_fasta_cmd = ['python ~/bin/JunkDrawer/fast_sweep2.py', filtered_fasta_file,
						clean2_fasta_file
						]
	cluster_ref_cmd = ['~/bin/cdhit/cd-hit -i', clean2_fasta_file, '-o', cluster_fasta_file,
						'-c 0.97 -d 0 -M 16000 -T 4'
						]
	final_ref_lineage_cmd = ['python ~/bin/JunkDrawer/build_ref_lineage.py',
								acc2taxid_lineage_file, cluster_fasta_file, run_path
								]
	clean3_fasta_cmd = ['python ~/bin/JunkDrawer/fast_sweep3.py', cluster_fasta_file,
						clean3_fasta_file, tax_id_file
						]
	build_ref_aln_cmd = ['muscle', '-in', clean3_fasta_file, '-out', final_aln_file]
	build_ref_hmm_cmd = ['hmmbuild --cpu 4', final_hmm_file, final_aln_file]
	build_ref_tree_cmd = ['FastTree -lg', final_aln_file, '>', final_tree_file]
	cp_tax_id_cmd = ['cp', tax_id_file, final_tax_id_file]

	#build_refpkg_cmd = ['~/bin/TreeSAPP/create_treesapp_ref_data.py -i', clean2_fasta_file,
	#					'-o', run_path,	'-c', run_name.rsplit('_', 1)[0],
	#					'-p 0.90 --cluster --trim_align -m prot -T 4 --headless'
	#					]
	if '.fasta' in file:
		cmds = [run_path, align_cmd, hmmbuild_cmd, hmmsearch_cmd, sto2fasta_cmd,
				clean_fasta_cmd, build_ref_lineage_cmd, sort_fasta_cmd,
				filter_fasta_cmd, build_refpkg_cmd]
	elif '.hmm' in file:
		cmds = [run_path, hmmsearch_cmd, sto2fasta_cmd,
				clean_fasta_cmd, build_ref_lineage_cmd, sort_fasta_cmd,
				filter_fasta_cmd, clean2_fasta_cmd, cluster_ref_cmd,
				final_ref_lineage_cmd, clean3_fasta_cmd, build_ref_aln_cmd,
				build_ref_hmm_cmd, build_ref_tree_cmd, cp_tax_id_cmd] # , build_refpkg_cmd]
	cmds_list.append(cmds)

pool = Pool(processes=int(num_jobs))                                                        
run_multi = pool.map(run_cmds, cmds_list)


