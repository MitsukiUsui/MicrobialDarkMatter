#!/usr/bin/env python3

import sys
import pathlib
import logging

import pandas as pd
from collections import defaultdict

ROOT_PATH = pathlib.Path().joinpath('../../').resolve()
sys.path.append(str(ROOT_PATH))
import mylib
from mylib.db import IDManager, get_connection, load_name2id
from mylib.df import read_mmseqs

LOGGER = logging.getLogger(__name__)
HID = IDManager("hits")


def main(mmseqs_fp, hits_fp, error_fp):
    LOGGER.info("Start extracting name2id from DB. Be patient...")
    con = get_connection()
    cds2id = load_name2id("cdss", default=-1, con=con)
    LOGGER.info("loaded {} cds names".format(len(cds2id)))
    refseq2id = load_name2id("refseqs", default=-1, con=con)
    LOGGER.info("loaded {} refseq names".format(len(refseq2id)))

    mmseqs_df = read_mmseqs(mmseqs_fp)
    LOGGER.info("loaded {} mmseqs hits".format(len(mmseqs_df)))
    mmseqs_df["hit_id"] = mmseqs_df.index.map(lambda x : HID.new())
    mmseqs_df["cds_id"] = mmseqs_df["qname"].map(cds2id)
    mmseqs_df["refseq_id"] = mmseqs_df["sname"].map(refseq2id)
    mmseqs_df["coverage"] = mmseqs_df["length"] / mmseqs_df["qlength"]

    msk = (mmseqs_df["cds_id"]!=-1) & (mmseqs_df["refseq_id"]!=-1)
    hits_df = mmseqs_df[msk][["hit_id", "cds_id", "refseq_id", "length", "identity", "coverage"]]
    hits_df.to_csv(hits_fp, index=False, sep='\t', header=None)
    LOGGER.info("saved {} records to {}".format(len(hits_df), hits_fp))

    error_df = mmseqs_df[~msk]
    error_df.to_csv(error_fp, index=False, sep='\t', header=None)
    LOGGER.info("saved {} records to {}".format(len(error_df), error_fp))


if __name__=="__main__":
    logging.basicConfig(level=logging.INFO, datefmt="%m/%d/%Y %I:%M:%S",
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    mmseqs_fp = "/nfs_share/mitsuki/MicrobialDarkMatter/search/result_3.m8.best"
    hits_fp = "./data/hits.tsv"
    error_fp = "./data/hits_error.tsv" # also output failed records for debuging
    main(mmseqs_fp, hits_fp, error_fp)
