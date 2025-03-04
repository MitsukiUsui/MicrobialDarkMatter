#!/usr/bin/env python3

import sys
import pathlib
import logging

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

ROOT_PATH = pathlib.Path().joinpath('../../').resolve()
sys.path.append(str(ROOT_PATH))
from mylib.db import get_session, Genome, Cds
from mylib.path import build_local_filepath

LOGGER = logging.getLogger(__name__)


def main(out_fp):
    session = get_session()
    genomes = session.query(Genome).all()
    LOGGER.info("found {} genomes to load".format(len(genomes)))

    records = []
    for genome in tqdm(genomes):
        cdss = session.query(Cds).filter_by(genome_id=genome.genome_id).all()
        cds2id = dict([(cds.cds_name, cds.cds_id) for cds in cdss])
        faa_fp = build_local_filepath(genome.genome_name, "faa")
        for seqrec in SeqIO.parse(faa_fp, "fasta"):
            cds_name = seqrec.id
            records.append({
                "cds_id": cds2id[cds_name],
                "seq": str(seqrec.seq),
            })
    out_df = pd.DataFrame(records, columns=["cds_id", "seq"])
    out_df.to_csv(out_fp, index=False, header=None, sep='\t')
    LOGGER.info("saved {} records to {}".format(len(out_df), out_fp))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, datefmt="%m/%d/%Y %I:%M:%S",
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    out_fp = "./data/sequences.tsv"
    main(out_fp)
