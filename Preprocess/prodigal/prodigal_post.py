#!/usr/bin/env python3

import sys
import pathlib
import logging

ROOT_PATH = pathlib.Path().joinpath('../../').resolve()
sys.path.append(str(ROOT_PATH))
import mylib
from mylib.gff import GffRecord

LOGGER = logging.getLogger(__name__)

"""
add cds_name attribute to Prodigal generated gff
"""

def main(genome_name, in_fp, out_fp):
    LOGGER.info("load {}".format(in_fp))
    cds_names = []
    lines = []
    with open(in_fp, "r") as f:
        for line in f:
            if line[0] == '#':
                lines.append(line.strip())
            else:
                record = GffRecord(line)
                assert record.type == "CDS"
                cds_name = "{}_{}".format(record.seqid, record.attributes["ID"].split('_')[1])
                record.attributes["cds_name"] = cds_name # add cds_name to attribute
                lines.append(str(record))
                cds_names.append(cds_name)
    assert len(set(cds_names)) == len(cds_names) # check if generated cds_names are unique

    with open(out_fp, 'w') as f:
        for line in lines:
            f.write("{}\n".format(line))
    LOGGER.info("output gff records with cds_name attribute to {}".format(out_fp))

if __name__=="__main__":
    logging.basicConfig(level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S",
                            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    genome_name = sys.argv[1]
    in_fp = sys.argv[2]
    out_fp = sys.argv[3]
    main(genome_name, in_fp, out_fp)

