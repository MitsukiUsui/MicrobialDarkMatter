#!/usr/bin/env python3

import sys
import pathlib
import logging

import pandas as pd

ROOT_PATH = pathlib.Path().joinpath('../../').resolve()
sys.path.append(str(ROOT_PATH))
import mylib
from mylib.path import build_clade_filepath, build_local_filepath
from mylib.db import get_connection

LOGGER = logging.getLogger(__name__)

def main(clade_name):
    con = get_connection()
    query = "SELECT genome_name FROM clades WHERE clade_name='{}'".format(clade_name)
    clades_df = pd.read_sql_query(query, con)
    LOGGER.info("found {} {} genomes".format(len(clades_df), clade_name))

    ortho_direc = pathlib.Path(build_clade_filepath(clade_name)).joinpath("./ortho/")
    ortho_direc.mkdir(parents=True, exist_ok=True)
    for genome_name in clades_df["genome_name"]:
        local_fp = pathlib.Path(build_local_filepath(genome_name, "faa"))
        input_fp = ortho_direc.joinpath(local_fp.name)
        if not(input_fp.exists()):
            input_fp.symlink_to(local_fp)
    LOGGER.info("prepared input directory. Please exec below:")
    cmd = "qsub ./orthofinder.sh {}".format(ortho_direc)
    print(cmd)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, datefmt="%m/%d/%Y %I:%M:%S",
                            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    clade_name = "Enterobacterales"
    main(clade_name)
