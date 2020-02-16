#!/usr/bin/env bash

OUT='.'
while getopts r:l:o option 
do
case "${option}" 
in
r) READS=${OPTARG};;
l) LENGTH=${OPTARG};;
o) OUT=${OPTARG};;
esac
done

case ${READS} in
    *.fastq.gz ) 
    BN=$(basename ${READS} .fastq.gz)
    zcat ${READS} | paste - - - - | awk -v L=${LENGTH} -F"\t" 'length($2)  >= L' | sed 's/\t/\n/g' | gzip > "${OUT}/${BN}_filtered.fastq.gz"
    ;;
    *.fastq)
    BN=$(basename ${READS} .fastq)
    cat ${READS} | paste - - - - | awk -v L=${LENGTH} -F"\t" 'length($2)  >= L' | sed 's/\t/\n/g' | gzip > "${BN}_filtered.fastq.gz"
    ;;
esac   
