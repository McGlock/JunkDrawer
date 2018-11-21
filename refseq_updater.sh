#!/usr/bin/env bash

# functions
function join_by { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }

# start script
workdir=/mnt/nfs/sharknado/LimsData/Hallam_Databases/raw/Protein/refseq

mkdir -p $workdir
db_list=( archaea bacteria fungi mitochondrion plasmid plastid protozoa viral )

for db in "${db_list[@]}"
do
	mkdir -p $workdir/$db
	db_dir=$workdir/$db
	wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/$db/*.protein.faa.gz -P $db_dir/
done

file_list=$( ls $workdir/*/*.protein.faa.gz)
DATE=`date +%Y-%m-%d`
final_file_name=$( join_by _ ${db_list[@]}).protein.$DATE.faa.gz
echo $final_file_name
cat $file_list > $workdir/$final_file_name
rm -fr $workdir/*/
