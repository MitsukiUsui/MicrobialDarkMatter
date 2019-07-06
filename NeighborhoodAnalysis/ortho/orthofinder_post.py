#!/usr/bin/env python3

import sys
import pathlib
import logging

import pandas as pd

ROOT_PATH = pathlib.Path().joinpath('../../').resolve()
sys.path.append(str(ROOT_PATH))
from mylib.path import build_clade_filepath

LOGGER = logging.getLogger(__name__)


def main(clade_name):
    ortho_direc = pathlib.Path(build_clade_filepath(clade_name)).joinpath("./ortho/")
    #    result_fp = list(ortho_direc.rglob("Orthogroups.csv"))[0]
    #    result_fp = "/home/mitsuki/afp/material/Enterobacterales/ortho/Results_Nov14/Orthogroups.csv"
    result_fp = "/home/mitsuki/afp/material/MGII_95/ortho/Results_Jan11/Orthogroups.csv"
    out_fp = ortho_direc.joinpath("{}.ortho".format(clade_name))

    result_df = pd.read_csv(result_fp, sep='\t', index_col=0, dtype=object)  # index: gene_name, column: genome_name
    LOGGER.info("loaded raw results from {}".format(result_fp))

    records = []
    for genome_name in result_df.columns:
        for gene_name, value in result_df[genome_name].dropna().iteritems():
            for cds_name in value.split(','):
                records.append({
                    "gene_name": gene_name,
                    "cds_name": cds_name.strip()
                })
    out_df = pd.DataFrame(records, columns=["gene_name", "cds_name"])
    out_df = out_df.sort_values(by=["gene_name", "cds_name"])
    out_df.to_csv(out_fp, index=False, sep='\t')
    LOGGER.info("saved parsed result to {}".format(out_fp))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, datefmt="%m/%d/%Y %I:%M:%S",
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    #    clade_name = "Enterobacterales"
    clade_name = "MGII"
    main(clade_name)
