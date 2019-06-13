#!/usr/bin/env python3

import sys
import pathlib
import logging

import pandas as pd
from collections import defaultdict

ROOT_PATH = pathlib.Path().joinpath('../../').resolve()
sys.path.append(str(ROOT_PATH))
import mylib
from mylib.db import IDManager, get_connection
from mylib.df import read_mmseqs

LOGGER = logging.getLogger(__name__)
HID = IDManager("hits")


def to_ddct(keys, vals, default):
    ddct = defaultdict(lambda: default)
    for key, val in zip(keys, vals):
        ddct[key] = val
    return ddct


def main(mmseqs_fp, hits_fp, error_fp):
    con = get_connection()
    query = "SELECT cds_id, cds_name FROM cdss;"
    cdss_df = pd.read_sql_query(query, con)
    cds2id = to_ddct(cdss_df["cds_name"], cdss_df["cds_id"], default=-1)
    LOGGER.info("load cdss")

    query = "SELECT refseq_id, refseq_name FROM refseqs;"
    refseqs_df = pd.read_sql_query(query, con)
    refseq2id = to_ddct(refseqs_df["refseq_name"], refseqs_df["refseq_id"], default=-1)
    LOGGER.info("load refseqs")

    mmseqs_df = read_mmseqs(mmseqs_fp)
    mmseqs_df["hit_id"] = mmseqs_df.index.map(lambda x : HID.new())
    mmseqs_df["cds_id"] = mmseqs_df["qname"].map(cds2id)
    mmseqs_df["refseq_id"] = mmseqs_df["sname"].map(refseq2id)
    mmseqs_df["coverage"] = mmseqs_df["length"] / mmseqs_df["length"]

    msk = (mmseqs_df["cds_id"]!=-1) & (mmseqs_df["refseq_id"]!=-1)
    hits_df = mmseqs_df[msk][["cds_id", "refseq_id", "length", "identity", "coverage"]]
    hits_df.to_csv(hits_fp, index=False, sep='\t', header=None)
    LOGGER.info("output {}".format(hits_fp))

    error_df = mmseqs_df[~msk]
    error_df.to_csv(error_fp, index=False, sep='\t', header=None)
    LOGGER.info("output {}".format(error_fp))


if __name__=="__main__":
    logging.basicConfig(level=logging.INFO, datefmt="%m/%d/%Y %I:%M:%S",
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    mmseqs_fp = "/home/mitsuki/data/mag/mmseqs/cdss30_REP_3.best.head"
    hits_fp = "./data/hits.tsv"
    error_fp = "./data/hits_error.tsv" # also output failed records for debuging
    main(mmseqs_fp, hits_fp, error_fp)
