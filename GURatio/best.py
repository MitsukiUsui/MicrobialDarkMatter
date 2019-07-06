#!/usr/bin/env python3

"""
Extract only the best MMSeqs hit for each query
"""

import sys
import pathlib

ROOT_PATH = pathlib.Path().joinpath('..').resolve()
sys.path.append(str(ROOT_PATH))
from mylib.df import read_mmseqs


def main(mmseqs_fp, best_fp):
    mmseqs_df = read_mmseqs(mmseqs_fp, dtype=object)
    best_df = mmseqs_df.drop_duplicates("qname")
    assert len(best_df) == len(set(mmseqs_df["qname"]))
    best_df.to_csv(best_fp, sep='\t', index=False, header=None)
    print("saved {}".format(best_fp))


if __name__ == "__main__":
    mmseqs_fp = "/nfs_share/mitsuki/MicrobialDarkMatter/search/result_3.m8"
    best_fp = mmseqs_fp + ".best"
    main(mmseqs_fp, best_fp)
