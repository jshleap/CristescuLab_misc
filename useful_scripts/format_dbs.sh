#!/usr/bin/env bash

# This script takes a fasta from NCBI or BOLD, and format it with taxonomic info

#______________________________________
# Copyright (C) 2018  Jose Sergio Hleap
#______________________________________

file=$1
typ=$2

add_lineages2accession2taxid(){
#$1 path to concatenation of accession2taxid
com='taxonkit lineage -i 3 {}| taxonkit reformat -i 4 > accession2taxid.lineages'
parallel $com  ::: <(cut -f1,2,3 $1)
}


process_NCBI(){
# Field 1 is the fasta file
# Field 2 is the path to the concatenation of accession2taxid
# Field 3 output filename
com="grep -w -Ff {1} {2} > $3"
grep '>' $1| sed 's/>//g'| cut -d' ' -f1 > accession.list
parallel $com ::: accession.list ::: $2
comm -13 <(cut -f 1 $3| sort -u) <(sort accession.list) > difference
parallel -a <(sort -u difference) ./get_failed_acc.sh > missing1.txt
}
#Get accessions from NCBI FASTA