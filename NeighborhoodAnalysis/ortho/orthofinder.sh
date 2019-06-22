#!/bin/bash
#$ -S /bin/bash
#$ -N orthofinder
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/ortho_$JOB_ID.out
#$ -e ./log/ortho_$JOB_ID.err
#$ -pe make 30

set -ue

direc=${1}
orthofinder -f ${direc} -S mmseqs -t 30
