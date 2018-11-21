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


# number of cores on system
num_jobs = sys.argv[4] # multiprocessing.cpu_count()
# set paths
home_path = sys.argv[1]
file_path = sys.argv[2]
refpkg_dir = sys.argv[3]
ref_db_path = sys.argv[5]
acc2taxid_lineage_file = sys.argv[6]

# create refpkg directory
refpkg_path = makedir(home_path, refpkg_dir)
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
for file in file_path:
	run_name = file.split('/')[-1].split('.')[0]
	run_path = makedir(refpkg_path, run_name)
	ref_aln_file = join(run_path, run_name + '.ref.aln')
	if '.fasta' in file:
		ref_hmm_file = join(run_path, run_name + '.ref.hmm')
	elif '.hmm' in file:
		ref_hmm_file = file
	recruit_aln_file = join(run_path, run_name + '.recruit.sto')
	recruit_fasta_file = join(run_path, run_name + '.recruit.fasta')
	clean_fasta_file = join(run_path, run_name + '.clean.fasta')
	exclude_acc_file = join(run_path, run_name + '.exclude.txt')
	filtered_fasta_file = join(run_path, run_name + '.filtered.fasta')
	align_cmd = ['muscle', '-in', file, '-out', ref_aln_file]
	hmmbuild_cmd = ['hmmbuild --cpu 4', ref_hmm_file, ref_aln_file]
	hmmsearch_cmd = ['hmmsearch --cpu 4 -A', recruit_aln_file, ref_hmm_file, ref_db_path]
	sto2fasta_cmd = ['seqmagick convert --ungap', recruit_aln_file,  recruit_fasta_file]
	clean_fasta_cmd = ['python ~/bin/fast_sweep.py', recruit_fasta_file, clean_fasta_file]
	build_ref_lineage_cmd = ['python ~/bin/build_ref_lineage.py', acc2taxid_lineage_file,
								clean_fasta_file, run_path
								]
	exclude_acc_cmd = ['seqmagick convert', clean_fasta_file, filtered_fasta_file,
						'--exclude-from-file', exclude_acc_file
						]
	build_refpkg_cmd = ['~/bin/TreeSAPP/create_treesapp_ref_data.py -i', filtered_fasta_file,
						'-o', run_path,	'-c', run_name.rsplit('_', 1)[0],
						'-p 0.90 --cluster --trim_align -m prot -T 4 --headless'
						]
	if '.fasta' in file:
		cmds = [run_path, align_cmd, hmmbuild_cmd, hmmsearch_cmd, sto2fasta_cmd,
				clean_fasta_cmd, build_ref_lineage_cmd, exclude_acc_cmd, build_refpkg_cmd]
	elif '.hmm' in file:
		cmds = [run_path, hmmsearch_cmd, sto2fasta_cmd,
				clean_fasta_cmd, build_ref_lineage_cmd, exclude_acc_cmd, build_refpkg_cmd]
	cmds_list.append(cmds)

#run_para = Parallel(n_jobs=int(num_jobs))(delayed(run_cmds)(cmds) for cmds in [cmds_list[0]])
pool = Pool(processes=int(num_jobs))                                                        
run_multi = pool.map(run_cmds, cmds_list)


