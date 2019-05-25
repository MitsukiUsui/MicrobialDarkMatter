#!/bin/bash

set -ue

genome_name=${1}
gz_fp=${2}
fna_fp=${3}
faa_fp=${4}
gff_fp=${5}

#if [ ! -e ${gff_fp} ]; then
    # generate .fna with parsable sequence_id
    gunzip -dc ${gz_fp} > ${gz_fp/.gz/}
    ./prodigal_pre.py ${genome_name} ${gz_fp/.gz/} ${fna_fp}
    rm ${gz_fp/.gz/}

    prodigal -i ${fna_fp} \
             -a ${faa_fp} \
             -o ${gff_fp} -f gff

    # assign cds_name to attributes in .gff
    ./prodigal_post.py ${genome_name} ${gff_fp} ${gff_fp}
#fi
