#!/usr/bin/env python3

"""
rename sequence_id in .fna
"""

import sys
import logging

from Bio import SeqIO

LOGGER = logging.getLogger(__name__)


def main(genome_name, in_fp, out_fp):
    LOGGER.info("load {}".format(in_fp))
    records = []
    for record in SeqIO.parse(in_fp, "fasta"):
        scaffold_name = "{}-{}".format(genome_name, record.id)  # rename to scaffold_id
        record.id = scaffold_name
        records.append(record)

    with open(out_fp, "w") as f:
        SeqIO.write(records, f, "fasta")
    LOGGER.info("saved renamed records to {}".format(out_fp))

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S",
                                format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    genome_name = sys.argv[1]
    in_fp = sys.argv[2]
    out_fp = sys.argv[3]
    main(genome_name, in_fp, out_fp)
