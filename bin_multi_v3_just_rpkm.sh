#!/bin/bash


if [ $# -lt 6 ]; then
        echo -e "\n\tUSAGE: assembly_list.txt reads_list.txt interleaved[y|n] output_dir prefix_name n_threads\n"
	echo "'assembly_list.txt' is a file listing the FASTA files containing contigs to be binned."
	echo "                    They will be concatenated, length-filtered and deduplicated prior to binning."
	echo "'reads_list.txt'    is a file listing the FASTQ files to be used for generating the abundance profiles."
	echo "                    Interleaved FASTQ will need a single file on one line; split FASTQ needs R1 and R2 tab-delimited."
	echo "                    A split-FASTQ file list can be created with '$ ls /path/to/*/trim/*pe.?.fq.gz | paste - - >reads_list.txt'"
	echo "'interleaved'       specifies whether the FASTQ files of each sample are split (n) or interleaved (y)."
	echo "'output_dir'        is the directory for writing output files."
	echo "'prefix_name'       is the prefix for naming output files. Must be less than 20 characters."
	echo "'n_threads'         is the number of threads to be used by bwa, samtools, and MetaBAT 2.11.2."
        exit
fi

input_contig=$1
reads_list=$2
fq_form=$3
prefix_dir=$4
sid=$5
# f_ext=$( echo $input_contig | awk -F. '{ print $NF }')
# sid=$( basename $input_contig .$f_ext)
n_thread=$6

min_length=1500
dedup_id=99
log_f=$prefix_dir/${sid}_gt${min_length}_BINME_log.txt
abunds=$prefix_dir/${sid}_gt${min_length}_abund_list.txt
bin_input_contigs=$prefix_dir/${sid}_gt${min_length}_dedup.fasta
output_dir=$prefix_dir/MetaBAT2_$sid/
bad_dir=$output_dir/bad_name_bins
rpkm_f=$prefix_dir/${sid}_binned.rpkm.csv
header_map=$output_dir/binned_header_map.tsv
checkm_stats=$output_dir/MetaBAT2_${sid}_min${min_length}_checkM_stdout.tsv

# Ensure the prefix is not too long - essential for compatibility with Prokka
if [ $(echo $sid| wc -c ) -gt 21 ]; then
	echo "ERROR: 'prefix_name' ($sid) must be shorter than 20 characters! Exiting now."
	exit
fi

# Ensure formatting of FASTQ files is correct
if [ $fq_form != 'y' ] && [ $fq_form != 'n' ]; then
	echo "ERROR: Expected a 'y' or 'n' for the interleaved argument. Cannot recognize '$fq_form'! "
	exit
fi
while read line
do
	nfq=$( echo "$line" | awk --field-separator="\t" '{ print NF }' )
	if [ $nfq -eq 2 ] && [ $fq_form == 'y' ]; then
		echo "ERROR: 2 FASTQs found when only one is expected in $reads_list."
		exit
	elif [ $nfq -eq 1 ] && [ $fq_form == 'n' ]; then
		echo "ERROR: 1 FASTQ found when 2 are expected in $reads_list."
		exit
	fi
done<$reads_list

printf "\n################################################## BIN_MULTI ##################################################\n"

printf "\nPARAMETERS:\n"
printf "\tPREFIX = $sid \n"
printf "\tOUTPUT DIRECTORY = $prefix_dir \n"
printf "\tMIN_CONTIG_LENGTH = $min_length \n"
printf "\tLOG = $log_f \n"
printf "\tTHREADS = $n_thread \n"
printf "\tMIN_DEDUP_ID = $dedup_id \n"
echo

if [ ! -d $prefix_dir ]; then
	mkdir $prefix_dir
fi

if [ -f $log_f ]; then
	rm $log_f
	touch $log_f
fi 

###################################################################################################
# Calculate RPKM for all binned contigs
###################################################################################################
printf "[STATUS] Creating BWA index for $output_dir/binned_sequences.fasta... "
cat $output_dir/MAGs/*.fa >$output_dir/binned_sequences.fasta
if [ ! -f $output_dir/binned_sequences.fasta ]; then
	echo "ERROR: Unable to concatenate bin FASTA files in $output_dir/"
	exit
fi
bwa index $output_dir/binned_sequences.fasta 1>$prefix_dir/bwa_stdout.txt 2>$prefix_dir/bwa_stderr.txt
printf "done.\n"

touch $rpkm_f
echo "Sample,Sequence,RPKM" >>$rpkm_f
while read line; do
	fq_sid=$( basename $line | sed 's/.fastq.*\|.fq.*//g' )
	printf "[STATUS] Aligning reads in $fq_sid ... "
	if [ $fq_form == 'y' ]; then
	bwa mem -t $n_thread -p $output_dir/binned_sequences.fasta $line 1>$prefix_dir/binned_sequences.$fq_sid.sam 2>$prefix_dir/bwa_stderr.txt
	else
		fwfq=$( echo "$line" | awk -F"\t" '{ print $1 }' )
		refq=$( echo "$line" | awk -F"\t" '{ print $2 }' )
	bwa mem -t $n_thread $output_dir/binned_sequences.fasta $fwfq $refq 1>$prefix_dir/binned_sequences.$fq_sid.sam 2>$prefix_dir/bwa_stderr.txt
	fi
	printf "done.\n"
	printf "[STATUS] Calculating RPKM... "
	rpkm -c $output_dir/binned_sequences.fasta -a $prefix_dir/binned_sequences.$fq_sid.sam -o $prefix_dir/$sid\_binned.$fq_sid.rpkm.csv 1>>$log_f 2>>$log_f
	cat $prefix_dir/$sid\_binned.$fq_sid.rpkm.csv | sed "s/^/$fq_sid,/g" >>$rpkm_f
	rm $prefix_dir/$sid\_binned.$fq_sid.rpkm.csv
	printf "done.\n"
done<$reads_list

rm $prefix_dir/binned_sequences.*.sam $output_dir/binned_sequences.fasta*

