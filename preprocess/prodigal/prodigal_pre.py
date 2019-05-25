#!/usr/bin/env python3

"""
rename sequence_id in .fna
"""

import sys
import logging

from Bio import SeqIO

logger = logging.getLogger(__name__)


def main(genome_name, original_fp, renamed_fp):
    logger.info("load {}".format(original_fp))
    with open(renamed_fp, "w") as f:
        for record in SeqIO.parse(original_fp, "fasta"):
            record.id = "{}-{}".format(genome_name, record.id)  #add  genome_id to head
            SeqIO.write(record, f, "fasta")
    logger.info("saved renamed records to {}".format(renamed_fp))

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S",
                                format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    genome_name = sys.argv[1]
    original_fp = sys.argv[2]
    renamed_fp = sys.argv[3]
    main(genome_name, original_fp, renamed_fp)
