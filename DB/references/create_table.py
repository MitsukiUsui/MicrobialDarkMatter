#!/usr/bin/env python3

import re
import sys
import pathlib
import logging

import pandas as pd

ROOT_PATH = pathlib.Path().joinpath('../../').resolve()
sys.path.append(str(ROOT_PATH))
from mylib.db import IDManager

LOGGER = logging.getLogger(__name__)
RID = IDManager("refseqs")
HEADER_PATTERN = re.compile(r"(.P_[0-9.]+) (MULTISPECIES:)?([^\[\]]*)(\[(.*)\])?")


def is_function_known(desc):
    desc = desc.lower()
    if "hypothetical protein" in desc or desc == "membrane protein":
        return False
    return True


def parse_header(header):
    record = dict()
    if isinstance(header, str):
        match = re.match(HEADER_PATTERN, header)
        if match:
            record["accession"] = match.group(1)
            record["description"] = match.group(3).strip()
            record["lca"] = match.group(5)
    return record


def main(stat_fp, out_fp):
    stat_df = pd.read_csv(stat_fp, sep="\t\t\t", header=None, names=["header", "length"], engine="python")
    LOGGER.info("loaded {} records from {}".format(len(stat_df), stat_fp))

    records = list(map(parse_header, stat_df["header"]))
    header_df = pd.DataFrame(records, columns=["accession", "description", "lca"])
    LOGGER.info("parsed {} header".format(len(header_df)))

    out_df = pd.DataFrame()
    out_df["refseq_id"] = stat_df.index.map(lambda x: RID.new())
    out_df["refseq_name"] = header_df["accession"]
    out_df["description"] = header_df["description"]
    out_df["lca"] = header_df["lca"]
    out_df["gk"] = 1
    out_df["fk"] = header_df["description"].map(lambda desc: is_function_known(desc)).astype(int)
    out_df.drop_duplicates("refseq_name", inplace=True)

    out_df.to_csv(out_fp, index=False, header=None, sep='\t')
    LOGGER.info("outputed {} records to {}".format(len(out_df), out_fp))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, datefmt="%m/%d/%Y %I:%M:%S",
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    stat_fp = "/home/mitsuki/GeneNeighborhoodAnalysis/GURatio/data/refseqs.stat"
    out_fp = "./data/refseqs.tsv"
    main(stat_fp, out_fp)
