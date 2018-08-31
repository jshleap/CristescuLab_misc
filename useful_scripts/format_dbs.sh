#!/usr/bin/env bash

# This code transforms a fasta file to Cristescu-lab blast format and to dada2 format
#______________________________________
# Copyright (C) 2018  Jose Sergio Hleap
#______________________________________

# Arguments:
#   1. Fasta file
#   2. Path to accession2taxid lineages
#   3. Number of available cores
#   4. BOLD or NCBI
# The accessiion2taxid lineages is a file mapping the accession numbers, tax id,and lineages

#get the fasta prefix
acc2tax=$2
cpus=$3

process_NCBI()
{
# 1. fasta
# 2. acc2tax
# 3. cpu
#get accession list
prefix=${1%%.fasta}
grep '>' $1| cut -d' ' -f 1| sed 's/>//' > ${prefix}.accessions
# get the lineages based on the accession numbers
LC_ALL=C awk 'FNR==NR{a[$1];next} ($1 in a)' ${prefix}.accesions $2 > ${prefix}.lineages
# get the mapping of accession number and stitle for renaming
format_them $1 ${prefix} $3
}

process_BOLD()
{
echo
}


format_them()
{
# 1. fasta
# 2. prefix
# 3. cpus
grep -a '>' $1 | sed 's/ /\t/'|sed 's/>//'| \
awk -F $'\t' ' {OFS = FS} { t = $1; $1 = $2; $2 = t; print; } ' > $2.keyvalue
cut -f2,5 $2.lineages |sed 's/$/;/'> $2.map
# rename the fasta file
seqkit -j $3 replace -p ' (.+)$' -r ' {kv}' -k $2.keyvalue $1 | \
seqkit -j $3 replace -p ' (.+)$' -r ' {kv}' -k $2.map > $2_lineages.fa
sed 's/[^ ]* />/' $2_lineages.fa > $2_lineages4dada.fa
}


if [ "$4" = 'NCBI' ]; then
    process_NCBI $1 $2 $3
else
    process_BOLD