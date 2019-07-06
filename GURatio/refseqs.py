#!/usr/bin/env python3

import sys
import pathlib

import pandas as pd

ROOT_PATH = pathlib.Path().joinpath('..').resolve()
sys.path.append(str(ROOT_PATH))
from mylib.ncbi import build_ftp_filepath
from mylib.path import build_local_filepath


def main(meta_fp, arg_fp):
    meta_df = pd.read_csv(meta_fp, sep='\t')
    records = []
    for _, row in meta_df.iterrows():
        fp = pathlib.Path(build_local_filepath(row["genome_name"])) / \
             pathlib.Path(build_ftp_filepath(row["ncbi_acc"], "faa")).name.replace(".gz", '')
        if not (fp.exists()):
            print("{} does not exist".format(fp))
            exit()
        records.append({"fp": fp})
    arg_df = pd.DataFrame(records)
    arg_df.to_csv(arg_fp, sep='\t', index=False, header=None)
    print("output {} records to {}".format(len(arg_df), arg_fp))


if __name__ == "__main__":
    meta_fp = "/home/mitsuki/GeneNeighborhoodAnalysis/dataCollection/data/meta/REP.tsv"
    arg_fp = "./arg/refseqs.arg"
    main(meta_fp, arg_fp)
