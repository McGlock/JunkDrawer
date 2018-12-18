#!/bin/bash

##
# This (when there are no options provided) is how the help statement is triggered and printed to stdout
# It doesn't check the order, every argument is positional - so don't mess it up!
# The rdp_data_path will need to be changed to a different directory depending on where these data live on your system
##
if [ $# -eq 0 ]; then
	echo -e "\n$0 itags_dir qiime_mapping_file.txt similarity run_name num_threads rdp_data_path"
	echo -e "\n\tThis script runs the qiime pipeline for iTags. It overwrites all previous runs."
	echo -e "itags_dir\n\tpath to directory containing all FASTQ files with illumina-generated rRNA sequences."
	echo -e "similarity\n\tThe sequence similarity threshold for pick_otus.py (e.g., 0.97)."
	echo -e "run_name\n\tPrefix to give files/directories unique names."
	echo -e "rdp_data_path\n\tPath to training and representative OTU set for RDP training (e.g.,"
	echo "/home/cmorganlang/bin/miniconda3/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus). [OPTIONAL]"
	echo -e "Order of operations:\n1. Check mapping file\n2. Demultiplex samples\n3. pick_otus.py"
	echo -e "4. pick_rep_set.py\n5. align_seqs.py using PyNAST\n6. identify_chimeric_seqs.py with UCHIME\n7. assign_taxonomy.py"
	exit
fi

# These are all explained in the help statement
itag_dir=$1
Fasting_map=$2
similarity=$3
prefix=$4
threads=$5
rdp_data=$6

##
# It gets a bit tricky here and a good estimate (and tight bounds) seem to help
# Georgetown sequences with the 515/806 primer set min = 250, max = 280
# Global soil atlas min = 410, max = 450
##
min_length=250
max_length=280

if [ -z $rdp_data ]; then
	rdp_data="/home/cmorganlang/bin/miniconda3/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/"
fi

# Remove the output directory if it already exists
if [ -d $prefix/ ]; then
	echo "WARNING: Output directory ($prefix) already exists and is being over-written."
	rm -r $prefix/
	mkdir $prefix/
else
	mkdir $prefix/
fi

# Check to ensure the similarity x holds for 0 < x < 1, otherwise exit
cond=$( echo "$similarity > 1" | bc)

if [ ! "$cond" ]; then
	echo "Bad similarity provided: $similarity"
	echo "similarity must be 0 < similarity < 1."
	exit
fi

# Check to make sure a binary called 'usearch61' is in your path
# This should really be USEARCH version 6.1 but we have been able to get away with some other versions. It is required by some QIIME steps
type usearch61 >/dev/null 2>&1 || { echo -e "ERROR:\nCould not find usearch61 executable in path. It is used for detecting chimeric reads"; echo "Please download it from http://www.drive5.com/usearch/download.html, rename the executable to usearch61, and add this to PATH"; exit; }

echo "Check mapping file"
if [ -d mapping_output ]; then
	rm -rf mapping_output
fi
validate_mapping_file.py -m $Fasting_map \
-o mapping_output \
--verbose  --disable_primer_check

cat $itag_dir/chunk*fastq > $itag_dir/merged.fastq
#cat $itag_dir/chunk*R1*fastq > $itag_dir/All_iTags_R1.fastq
#cat $itag_dir/chunk*R2*fastq > $itag_dir/All_iTags_R2.fastq

##
# The fastq_minovlen and minhsp must be equal otherwise it will take the higher of the two
# Feel free to change those parameters, though it likely won't change much within a few basepairs
##
# echo "Merging forward and reverse read pairs"
# usearch -fastq_mergepairs $itag_dir/All_iTags_R1.fastq -reverse $itag_dir/All_iTags_R2.fastq \
# -alnout aln.txt \
# -fastqout $itag_dir/usearch_merged.fastq \
# -report $prefix\_fastq_mergepairs_report.txt \
# -tabbedout $prefix\_fastq_mergepairs_table.txt \
# -threads $threads \
# -fastq_maxdiffs 10 \
# -fastq_minovlen 12 \
# -minhsp 12 
# exit
# echo "Compressing concatenated FASTQ files"
# rm $itag_dir/All_iTags_R?.fastq

echo "Extracting barcodes"
extract_barcodes.py -f $itag_dir/merged.fastq -l 14

echo "Demultiplexing. Sequences with new headers will be in $prefix\_split_library_output/seqs.fna "
if [ -d $prefix\_split_library_output ]; then
	rm -rf $prefix\_split_library_output
fi
split_libraries_fastq.py -m $Fasting_map \
-i reads.fastq \
-b barcodes.fastq \
-o $prefix\_split_library_output \
--max_barcode_errors=0 \
--barcode_type=6 \
--min_per_read_length_fraction=0.5 \
--sequence_max_n=10 \
--phred_offset=33 \
-v

num_raw_seqs=$( grep -c "^>" $prefix\_split_library_output/seqs.fna )
echo "Number of sequences in demultiplexed file = $num_raw_seqs"

echo "Checking for primers"
usearch -search_oligodb $prefix\_split_library_output/seqs.fna \
-db /mnt/nfs/sharknado/LimsData/Hallam_Databases/primer_db.fa \
-strand both \
-threads $threads \
-userout primer_search.txt \
-userfields query+target+qstrand+diffs+tlo+thi+trowdots

num_primer_positive=$(cat primer_search.txt | gawk '{ print $1 }' | sort -u | wc -l)
echo "Number of reads with primers attached = $num_primer_positive"
if [ $( echo "scale=2; $num_primer_positive/$num_raw_seqs > 0.90" | bc ) == 1 ]; then
	echo "Trimming primers since >90% of sequences contain them"
	usearch -fastx_truncate $prefix\_split_library_output/seqs.fna -stripleft 19 -stripright 20 -fastaout $prefix\_split_library_output/seqs_stripped.fna
	mv $prefix\_split_library_output/seqs_stripped.fna $prefix\_split_library_output/seqs.fna
else
	echo "Skipping primer trimming!"
fi

echo "Sorting the demultiplexed sequences by length"
echo "Fragments shorter than $min_length or longer than $max_length will be removed"
usearch -sortbylength $prefix\_split_library_output/seqs.fna -fastaout $prefix\_split_library_output/seqs_sorted.fasta \
-maxseqlength $max_length -minseqlength $min_length

echo "Identifying chimeras using UCHIME"
if [ -d $prefix\_chimeric_seqs ]; then
	rm -rf $prefix\_chimeric_seqs
fi
mkdir $prefix\_chimeric_seqs/

usearch -uchime_ref $prefix\_split_library_output/seqs_sorted.fasta \
-db /mnt/nfs/sharknado/LimsData/Hallam_Databases/rdp_gold.fa \
-uchimeout $prefix\_chimeric_seqs/results.uchime \
-chimeras $prefix\_chimeric_seqs/chimeras.fna \
-nonchimeras $prefix\_split_library_output/seqs_chimeras_filtered.fna \
-threads $threads \
-strand plus

grep "^>" $prefix\_chimeric_seqs/chimeras.fna | sed 's/^>//g' >$prefix\_chimeric_seqs/chimeras.txt

if [ -d $prefix\_picked_otus ]; then
	rm -rf $prefix\_picked_otus
fi

echo "Picking OTUs de novo with sumaclust"
pick_otus.py -i $prefix\_split_library_output/seqs_chimeras_filtered.fna \
-o $prefix\_picked_otus \
-v \
--otu_picking_method=sumaclust \
--sumaclust_exact \
--trie_prefilter \
--enable_rev_strand_match \
--threads $threads \
--similarity=$similarity

#echo "Picking OTUs via closed-reference for PiCRUSt"
#echo "pick_otus:enable_rev_strand_match True"  >> $prefix/otu_picking_params_97.txt
#echo "pick_otus:similarity 0.97" >> $prefix/otu_picking_params_97.txt
#pick_closed_reference_otus.py -i $prefix\_split_library_output/seqs_chimeras_filtered.fna \
#-o $prefix\_closed_ref_picked_otus \
#-p $prefix/otu_picking_params_97.txt \
#-r $rdp_data/rep_set/97_otus.fasta \
#-t $rdp_data/taxonomy/97_otu_taxonomy.txt \
#--parallel --jobs_to_start=$threads

echo "Picking representative set of sequences"
if [ -f $prefix\_repSet.fna ]; then
	rm $prefix\_repSet.fna
fi
pick_rep_set.py -i $prefix\_picked_otus/seqs_chimeras_filtered_otus.txt \
-f $prefix\_split_library_output/seqs_chimeras_filtered.fna \
-o $prefix\_repSet.fna

echo "Assign taxonomy with rdp classifier"
if [ -d $prefix\_assigned_taxa ]; then
	rm -rf $prefix\_assigned_taxa
fi
parallel_assign_taxonomy_rdp.py -i $prefix\_repSet.fna \
--rdp_max_memory=8000 \
-c 0.85 \
-t $rdp_data/taxonomy/97_otu_taxonomy.txt \
-r $rdp_data/rep_set/97_otus.fasta \
-o $prefix\_assigned_taxa \
--jobs_to_start=$threads

echo "Making OTU table"
if [ -f $prefix\_otu_table.biom ]; then
	rm $prefix\_otu_table.biom
fi
sid=$( basename $prefix )
make_otu_table.py \
--taxonomy $prefix\_assigned_taxa/$sid\_repSet_tax_assignments.txt \
-i $prefix\_picked_otus/seqs_chimeras_filtered_otus.txt \
-o $prefix\_otu_table.biom \
-e $prefix\_chimeric_seqs/chimeras.txt

echo "Converting otu_table.biom to Tab-separated format"
if [ -f $prefix\_otu_table.tsv ]; then
        rm $prefix\_otu_table.tsv
fi
biom convert -i $prefix\_otu_table.biom \
-o $prefix\_otu_table.tsv \
--to-tsv \
--table-type="OTU table"

map_taxonomy_to_OTUs.py -t $prefix\_assigned_taxa/$sid\_repSet_tax_assignments.txt \
-p $prefix\_otu_table.tsv \
-o $prefix\_otu_table_TaxaAssigned.tsv

echo "Filtering singletons from OTU table"
if [ -f $prefix\_otu_table_noSingletons.biom ]; then
        rm $prefix\_otu_table_noSingletons.biom
fi
filter_otus_from_otu_table.py -i $prefix\_otu_table.biom \
-o $prefix\_otu_table_noSingletons.biom \
--min_count 2 \
--min_samples 1

echo "Converting $prefix\_otu_table_noSingletons.biom to Tab-separated format"
if [ -f $prefix\_otu_table_noSingletons.tsv ]; then
        rm $prefix\_otu_table_noSingletons.tsv
fi
biom convert -i $prefix\_otu_table_noSingletons.biom \
-o $prefix\_otu_table_noSingletons.tsv \
--to-tsv \
--table-type="OTU table"

map_taxonomy_to_OTUs.py -t $prefix\_assigned_taxa/$sid\_repSet_tax_assignments.txt \
-p $prefix\_otu_table_noSingletons.tsv \
-o $prefix\_otu_table_noSingletons_TaxaAssigned.tsv

mv $prefix\_* $prefix/
rm barcodes.fastq reads.fastq
rm -r mapping_output/
mv primer_search.txt $prefix/

echo "Finished QIIME pipeline for $prefix"
