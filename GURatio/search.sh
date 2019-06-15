#!/bin/bash
#$ -S /bin/bash
#$ -N search
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/search_$JOB_ID.out
#$ -e ./log/search_$JOB_ID.err
#$ -l mem_free=100G
#$ -pe make 8

set -ue

threads=8
query_fp="./data/queries.faa"
target_fp="./data/refseqs.faa"
result_fp="/nfs_share/mitsuki/GeneNeighborhoodAnalysis/search/result.m8.head"

mmseqs_direc=/dev/shm/${USER}/mmseqs/`date +%s`
tmp_direc=${mmseqs_direc}/tmp
mkdir -p ${mmseqs_direc}
mkdir -p ${tmp_direc}
query_db=${mmseqs_direc}/queryDB
target_db=${mmseqs_direc}/targetDB
result_db=${mmseqs_direc}/resultDB

mmseqs createdb ${query_fp} ${query_db}
mmseqs createdb ${target_fp} ${target_db}
mmseqs search ${query_db} ${target_db} ${result_db} ${tmp_direc} --threads ${threads} -s 1 --max-seqs 20
mmseqs convertalis ${query_db} ${target_db} ${result_db} ${result_fp} --threads ${threads}
