#!/bin/bash

#NOTE: script assumes amos has been installed in the ~/bin directory

if [[ $# -eq 0 ]] ; then
    echo "./run_minimus fasta_list.txt prefix_string"
    echo "NOTE: by default this script is looking for high-confidence overlaps of 200bp @ 95%ID."
    echo "The intermediate files are deleted by default and the unscaffolded contigs are appended to
the minimus2 scaffolds (found in prefix_minimus2.fasta)"
    exit
fi

fasta_list=$1
prefix=$2

if [ -f $prefix\_unmerged.fasta ]; then
	echo "ERROR: $prefix\_unmerged.fasta already exists. Not overwriting!"
	exit
else
	touch $prefix\_unmerged.fasta
fi

numContigs=0

#usage: concatenate_FAs.py [-h] -o OUTPUT [-l FASTA_LIST]
#
#Script to concatenate multiple FASTA files and ensure no redundantheaders
#exist by prepending file names on the header.
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -o OUTPUT, --output OUTPUT
#                        The output FASTA file
#  -l FASTA_LIST, --fasta_list FASTA_LIST
#                        A list of FASTA files to concatenate with paths
#                        included


while read line
do
	if [ $numContigs -eq 0 ]; then
		numContigs=$( grep -c "^>" $line)
	fi
	if [ ! -f $line ]; then
		echo "ERROR: file $line does not exist!"
		exit
	fi
done<$fasta_list

echo "Concatenating assemblies into $prefix\_unmerged.fasta"
concatenate_FAs.py -o $prefix\_unmerged.fasta -l $fasta_list

total_contigs=$( grep -c "^>" $prefix\_unmerged.fasta)
echo "Reference assembly contains $numContigs / $total_contigs"

echo "Building AMOS database"
~/bin/amos-3.1.0/bin/toAmos -s $prefix\_unmerged.fasta -o $prefix.afg

echo "Running minimus2"
~/bin/amos-3.1.0/bin/minimus2 $prefix -D REFCOUNT=$numContigs -D OVERLAP=200 -D MINID=95
cat $prefix.fasta $prefix.singletons.seq >$prefix\_minimus2.fasta
rm -r $prefix.runAmos.log $prefix.afg $prefix.OVL $prefix.singletons \
$prefix.contig $prefix.ovl $prefix.singletons.seq $prefix.coords $prefix.qry.seq \
$prefix.delta $prefix.bnk $prefix.ref.seq $prefix.fasta
