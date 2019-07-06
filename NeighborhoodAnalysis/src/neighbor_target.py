#!/usr/bin/env python3

"""
Given a set of neighborhood relationships, calculate neighborhood scores for the given scoring method.
This script were used to calculate RMSE for each scoring method
"""

import sys
import pathlib
import logging
import argparse

import pandas as pd

ROOT_PATH = pathlib.Path().joinpath('../../').resolve()
sys.path.append(str(ROOT_PATH))
from mylib.db import CdsDAO, load_genome_names_by_clade_name, load_cdss_by_genome_names
from mylib.path import build_clade_filepath
from .neighborlib import NeighborhoodMatrix, set_gene_name, set_split
from .scorelib import score_naive, score_independent, score_conditional

LOGGER = logging.getLogger(__name__)


def detect_edges_target(origin_gene_name, neighbor_gene_names, score_method, cdsDAO):
    neighbor_matrix = NeighborhoodMatrix(origin_gene_name, cdsDAO)
    records = []
    for neighbor_gene_name in neighbor_gene_names:
        indicator_matrix = neighbor_matrix.to_indicator_matrix(neighbor_gene_name)
        if score_method == "naive":
            score = score_naive(indicator_matrix)
        elif score_method == "independent":
            score = score_independent(indicator_matrix)
        elif score_method == "conditional":
            score = score_conditional(indicator_matrix)
        records.append({
            "x": origin_gene_name,
            "y": neighbor_gene_name,
            "score": score
        })
    return records


def main(args):
    genome_names = load_genome_names_by_clade_name(args.clade_name)
    LOGGER.info("loaded {} {} genomes".format(len(genome_names), args.clade_name))
    cdss = load_cdss_by_genome_names(genome_names)
    LOGGER.info("loaded {} cdss".format(len(cdss)))

    ortho_fp = pathlib.Path(build_clade_filepath(args.clade_name)).joinpath("./ortho/{}.ortho".format(args.clade_name))
    cdss = set_gene_name(cdss, ortho_fp)
    LOGGER.info("loaded orthology from {}".format(ortho_fp))

    if args.split_fp:
        cdss = set_split(cdss, args.split_fp)
        LOGGER.info("loaded simulated segmentation from {}".format(args.split_fp))

    neighbor_df = pd.read_csv(args.neighbor_fp, comment='#')
    LOGGER.info("loaded {} relationships from {}".format(len(neighbor_df), args.neighbor_fp))

    records = []
    cdsDAO = CdsDAO(cdss)
    for origin_gene_name, sub_df in neighbor_df.groupby("x"):
        neighbor_gene_names = sub_df["y"]
        records += detect_edges_target(origin_gene_name, neighbor_gene_names, args.score_method, cdsDAO)

    out_df = pd.DataFrame(records, columns=["x", "y", "score"])
    out_df.to_csv(args.out_fp, sep='\t', index=False)
    LOGGER.info("saved results to {}".format(args.out_fp))


if __name__=="__main__":
    logging.basicConfig(level=logging.INFO, datefmt="%m/%d/%Y %I:%M:%S",
                            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--clade_name", required=True, help="clade_name")
    parser.add_argument("--score_method", required=True, choices=["naive", "independent", "conditional"])
    parser.add_argument("--out_fp", required=True)
    parser.add_argument("--neighbor_fp", required=True, help=".neighbor to follow")
    parser.add_argument("--split_fp", help="simulated segment map")
    args = parser.parse_args()
    main(args)
