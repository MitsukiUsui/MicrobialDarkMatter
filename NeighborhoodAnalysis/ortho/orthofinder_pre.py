#!/usr/bin/env python3

import sys
import pathlib
import logging

import pandas as pd

ROOT_PATH = pathlib.Path().joinpath('../../').resolve()
sys.path.append(str(ROOT_PATH))
from mylib.path import build_clade_filepath, build_local_filepath
from mylib.db import load_genome_names_by_clade_name

LOGGER = logging.getLogger(__name__)

def main(clade_name):
    genome_names = load_genome_names_by_clade_name(clade_name)
    LOGGER.info("found {} {} genomes".format(len(genome_names), clade_name))

    ortho_direc = pathlib.Path(build_clade_filepath(clade_name)).joinpath("./ortho/")
    ortho_direc.mkdir(parents=True, exist_ok=True)
    LOGGER.info("created input directory at {}".format(ortho_direc))

    for genome_name in genome_names:
        local_fp = pathlib.Path(build_local_filepath(genome_name, "faa"))
        input_fp = ortho_direc.joinpath(local_fp.name)
        if not(input_fp.exists()):
            input_fp.symlink_to(local_fp)
    LOGGER.info("prepared input directory")
    LOGGER.info("please exec below:")
    cmd = "qsub ./orthofinder.sh {}".format(ortho_direc)
    print(cmd)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, datefmt="%m/%d/%Y %I:%M:%S",
                            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    clade_name = "Enterobacterales"
    main(clade_name)
