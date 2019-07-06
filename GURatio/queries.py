#!/usr/bin/env python3

import sys
import pathlib

import pandas as pd

ROOT_PATH = pathlib.Path().joinpath('..').resolve()
sys.path.append(str(ROOT_PATH))
from mylib.path import build_local_filepath
from mylib.db import get_session, Genome, Cds


def main(arg_fp):
    session = get_session()
    genomes = session.query(Genome).all()
    print("found {} genomes".format(len(genomes)))

    records = []
    for genome in genomes:
        fp = pathlib.Path(build_local_filepath(genome.genome_name, "faa"))
        if not (fp.exists()):
            print("{} does not exist".format(fp))
            exit()
        records.append({"fp": fp})
    arg_df = pd.DataFrame(records)
    arg_df.to_csv(arg_fp, sep='\t', index=False, header=None)
    print("saved {} records to {}".format(len(arg_df), arg_fp))


if __name__ == "__main__":
    arg_fp = "./arg/queries.arg"
    main(arg_fp)
