#!/usr/bin/env python3

"""
Detect all gene neighborhood relationships under given thresholds
"""

import argparse
import logging
import pathlib
import sys
from collections import Counter

import pandas as pd
from ete3 import PhyloTree

ROOT_PATH = pathlib.Path().joinpath('../../').resolve()
sys.path.append(str(ROOT_PATH))
from mylib.db import CdsDAO, load_genome_names_by_clade_name, load_cdss_by_genome_names
from mylib.path import build_clade_filepath
from neighborlib import NeighborhoodMatrix, set_gene_name_to_cdss, set_split_to_cdss, calc_bls
from scorelib import score_naive, score_independent, score_conditional

LOGGER = logging.getLogger(__name__)
THRESH = {
    "SIZE": 10,  # lower limit for gene size
    "SCORE": 0.8  # lower limit for neighborhood score to be reported
}


def find_most_common_position(positions):
    def get_relationship(offset, is_forward):
        if is_forward:
            return "cooriented"
        elif offset < 0:
            return "divergent"
        else:
            return "convergent"

    assert len(positions) > 0
    c = Counter([(pos.offset, pos.is_forward) for pos in positions])
    (offset, is_forward), freq = c.most_common()[0]
    relationship = get_relationship(offset, is_forward)

    return {
        "top_offset": offset,
        "top_relationship": relationship,
        "top_ratio": freq / len(positions)
    }


def detect_edges_all(origin_gene_name, score_method, cdsDAO, tree=None):
    neighbor_matrix = NeighborhoodMatrix(origin_gene_name, cdsDAO)

    records = []
    if neighbor_matrix.shape[0] < THRESH["SIZE"]:
        LOGGER.debug("too small ortholog size = {}".format(neighbor_matrix.shape[0]))
        return records

    neighbor_gene_names = sorted(neighbor_matrix.get_neighbor_gene_names())
    LOGGER.debug("found {} candidate neighbor genes".format(len(neighbor_gene_names)))
    for neighbor_gene_name in neighbor_gene_names:
        indicator_matrix = neighbor_matrix.to_indicator_matrix(neighbor_gene_name)
        if score_method == "naive":
            score = score_naive(indicator_matrix)
        elif score_method == "independent":
            score = score_independent(indicator_matrix)
        elif score_method == "conditional":
            score = score_conditional(indicator_matrix)

        if score >= THRESH["SCORE"]:
            positions = neighbor_matrix.get_positions_by_gene_name(neighbor_gene_name)
            found_origin_names = set(map(lambda pos: pos.origin_name, positions))
            found_genome_names = set(map(lambda pos: pos.origin_name.split('-')[0], positions))
            record = {
                "x": origin_gene_name,
                "y": neighbor_gene_name,
                "score": score,
                "score_naive": score_naive(indicator_matrix),
                "total": neighbor_matrix.shape[0],
                "found": len(found_origin_names),
                "bls": calc_bls(found_genome_names, tree) if tree else -1
            }
            record.update(find_most_common_position(positions))
            records.append(record)
    LOGGER.info("found {} records".format(len(records)))
    return records


def main(args):
    genome_names = load_genome_names_by_clade_name(args.clade_name)
    LOGGER.info("loaded {} {} genomes".format(len(genome_names), args.clade_name))
    cdss = load_cdss_by_genome_names(genome_names)
    LOGGER.info("loaded {} cdss".format(len(cdss)))

    ortho_fp = pathlib.Path(build_clade_filepath(args.clade_name)).joinpath("./ortho/{}.ortho".format(args.clade_name))
    ortho_df = pd.read_csv(ortho_fp, sep='\t')
    cdss = set_gene_name_to_cdss(cdss, ortho_fp)
    LOGGER.info("loaded orthology from {}".format(ortho_fp))

    if args.split_fp:
        cdss = set_split_to_cdss(cdss, args.split_fp)
        LOGGER.info("loaded simulated segmentation from {}".format(args.split_fp))

    tree = None
    if args.tree_fp:
        tree = PhyloTree(args.tree_fp, format=1)
        LOGGER.info("loaded phylogenetic tree from {}".format(args.tree_fp))

    records = []
    cdsDAO = CdsDAO(cdss)
    gene_names = sorted(set(ortho_df["gene_name"]))
    #    gene_names = list(gene_names)[:100]
    LOGGER.info("found {} genes to search".format(len(gene_names)))
    for origin_gene_name in gene_names:
        LOGGER.info("start {}".format(origin_gene_name))
        records += detect_edges_all(origin_gene_name, args.score_method, cdsDAO, tree)

    out_df = pd.DataFrame(records, columns=["x", "y", "score", "score_naive", "total", "found", "bls",
                                            "top_offset", "top_relationship", "top_ratio"])
    out_df.to_csv(args.out_fp, sep='\t', index=False)
    LOGGER.info("saved results to {}".format(args.out_fp))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, datefmt="%m/%d/%Y %I:%M:%S",
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--clade_name", required=True, help="clade_name")
    parser.add_argument("--score_method", required=True, choices=["naive", "independent", "conditional"])
    parser.add_argument("--out_fp", required=True)
    parser.add_argument("--split_fp", help="simulated segment map")
    parser.add_argument("--tree_fp", help="newick tree")
    args = parser.parse_args()
    main(args)
