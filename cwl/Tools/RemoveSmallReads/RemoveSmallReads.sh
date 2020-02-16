#!/usr/bin/env bash

READS=$1
LENGTH=500
OUT='.'

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
