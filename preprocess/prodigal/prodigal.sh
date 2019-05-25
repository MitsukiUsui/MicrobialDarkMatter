#!/bin/bash

set -ue

genome_name=${1}
gz_fp=${2}
fna_fp=${3}
faa_fp=${4}
gff_fp=${5}

#if [ ! -e ${gff_fp} ]; then
    # generate .fna with parsable scaffold_name
    gunzip -dc ${gz_fp} > ${fna_fp}
    ./prodigal_pre.py ${genome_name} ${fna_fp} ${fna_fp}

    prodigal -i ${fna_fp} \
             -a ${faa_fp} \
             -o ${gff_fp} -f gff

    # assign cds_name to attributes in .gff
    ./prodigal_post.py ${genome_name} ${gff_fp} ${gff_fp}
#fi
