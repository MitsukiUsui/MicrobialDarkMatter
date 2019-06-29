#!/usr/bin/env python3

import sys
import pathlib
import logging
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd

ROOT_PATH = pathlib.Path().joinpath('../../').resolve()
sys.path.append(str(ROOT_PATH))
from mylib.db import CdsDAO, load_genome_names_by_clade_name, load_cdss_by_genome_names
from mylib.path import build_clade_filepath
from neighborlib import scan_neighborhoods, set_gene_name, set_split, output_neighbor_df
from scoring import score_naive, score_independent, score_conditional

LOGGER = logging.getLogger(__name__)


def find_most_common_position(positions):
    def get_relationship(offset, direc):
        if direc == True:
            return "cooriented"
        elif offset < 0:
            return "divergent"
        else:
            return "convergent"

    assert len(positions) > 0
    c = Counter([(pos.offset, pos.direc) for pos in positions])
    (offset, direc), freq = c.most_common()[0]
    relationship = get_relationship(offset, direc)

    dct = {
        "top_offset": offset,
        "top_relationship": relationship,
        "top_ratio": freq / len(positions)
    }
    return dct


def detect_edges_all(origin_gene_name, score_method, cdsDAO, tree=None):
	obs_mat, gene2positions = scan_neighborhoods(origin_gene_name, cdsDAO)
	records = []
    for neighbor_gene_name, positions in sorted(gene2positions.items(), key=lambda x:x[0]):
        if neighbor_gene_name == origin_gene_name:
            continue

        score_naive = scoring.score_naive(obs_mat, positions)
        if score_method == "naive":
            score = scoring.score_naive(obs_mat, positions)
        elif score_method == "independent":
            score = scoring.score_independent(obs_mat, positions)
        elif score_method == "conditional":
            score = scoring.score_conditional(obs_mat, positions)

        if score >= THRESH["SCORE"]:
            found_cds_ids = set([origin_cds_ids[pos.i] for pos in positions])
            genome_names = set([cdss[cds_id].cds_name.split('-')[0] for cds_id in found_cds_ids])
            bls = calc_bls(genome_names, tree) if tree else -1

            records.append({
                "x": origin_gene_name,
                "y": neighbor_gene_name,
                "total": len(origin_cds_ids),
                "found": len(found_cds_ids),
                "score": score,
                "score_naive": score_naive,
                "bls": bls,
            })
            dct.update(find_most_common_position(positions))
            dct_lst.append(dct)
    return records


def main(args):
    LOGGER.info("Gene Neighborhood with Missing data (GNM) Scoring module")

    #load genomes & cdss from SQLite Database
    genome_names = [line.strip() for line in open(args.genome_fp, 'r')]
    genomes = load_genomes(genome_names)
    LOGGER.info("loaded {} genomes".format(len(genomes)))
    cdss = load_cdss(genomes, to_dict=True)
    LOGGER.info("loaded {} cdss".format(len(cdss)))

    # load and set orthology info
    ortho_df = pd.read_csv(args.ortho_fp)
    cds2id = dict([(cds.cds_name, cds.cds_id) for cds in cdss.values()]) #key: cds_name, val: cds_id
    ortho_df["cds_id"] = ortho_df["cds_name"].map(lambda x:cds2id[x])
    gene2ids = defaultdict(list) #key: gene_name, val: list of cds_ids
    for gene_name, cds_id in zip(ortho_df["gene_name"], ortho_df["cds_id"]):
        gene2ids[gene_name].append(cds_id)
        cdss[cds_id].gene_name = gene_name #remember some cdss lack .gene_name attribute
    LOGGER.info("loaded orthology from {}".format(args.ortho_fp))

    #read tree for bls calculation
    if args.tree_fp:
        tree = PhyloTree(args.tree_fp, format=1)
        LOGGER.info("loaded phylogenetic tree from {}".format(args.tree_fp))
    else:
        tree = None

    #split
    if args.split_fp:
        cdss = set_split(cdss, args.split_fp)
        LOGGER.info("loaded simulated segmentation from {}".format(args.split_fp))

    #search edges one by one
    edges = []
#    gene_names = sorted(set(ortho_df["gene_name"]))
    gene_names = ["OG{:07}".format(i) for i in (1105, 1130, 1161, 1162, 1211, 1226, 994)]
    LOGGER.info("found {} genes to search".format(len(gene_names)))
    for gene_name in gene_names:
        LOGGER.info("start {}".format(gene_name))
        edges_ = search_edges(args.score_method, gene_name, gene2ids[gene_name], cdss, tree)
        LOGGER.info("found {} edges".format(len(edges_)))
        edges += edges_

    #output edges as DataFrame with comments
    edge_df = pd.DataFrame(edges)
    edge_df = edge_df[["x", "y", "total", "found", "score", "score_naive", "bls", "top_offset", "top_relationship", "top_ratio"]]
    output_edge_df(edge_df, args.out_fp)
    LOGGER.info("DONE: output {}".format(args.out_fp))


if __name__=="__main__":
	logging.basicConfig(level=logging.INFO, datefmt="%m/%d/%Y %I:%M:%S",
                            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--clade_name", required=True, help="clade_name")
    parser.add_argument("--score_method", required=True, choices=["naive", "independent", "conditional"])
    parser.add_argument("--neighbor_fp", required=True, help=".neighbor to follow")
    parser.add_argument("--out_fp", required=True)
    parser.add_argument("--split_fp", help="simulated segment map")
    parser.add_argument("--tree_fp", help="newick tree")
    args = parser.parse_args()
    main(args)
