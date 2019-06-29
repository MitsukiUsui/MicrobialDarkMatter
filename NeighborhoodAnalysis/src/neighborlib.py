#!/usr/bin/env python3

from logging import getLogger
from collections import defaultdict, Counter

import numpy as np
import pandas as pd

LOGGER = getLogger(__name__)


class MatrixPosition:
    """
    store position in NeighborhoodMatrix (NxM)
    """

    def __init__(self, i, j, offset, direction, cds_name, gene_name):
        self.i = i
        self.j = j
        self.offset = offset
        self.direction = direction
        self.cds_name = cds_name
        self.gene_name = gene_name

    def __lt__(self, other):
        assert self.i == other.i
        return abs(self.offset) < abs(other.offset)

    def __repr__(self):
        return "<MatrixPosition({},{})>".format(self.i, self.j)


class NeighborhoodMatrix:
    DIST = 5

    def __init__(self, origin_gene_name, cdsDAO):
        """
        initialize the following data structures:

        self.matrix: 2D matrix of MatrixPosition (HxW), where H = #origin-cdss and W = 2 * DIST + 1
        self.gene2positions: key: gene_name, val: list of MatrixPosition
        """

        # initialize matrix & gene2positions
        matrix = []
        gene2positions = defaultdict(list)
        for origin_cds in cdsDAO.get_cdss_by_gene_name(origin_gene_name):
            row = []
            for offset in range(-self.DIST, self.DIST+1):
                neighbor_cds = cdsDAO.get_neighbor_cds(origin_cds, offset)
                if neighbor_cds is None:
                    pos = None
                else:
                    pos = MatrixPosition(len(matrix), len(row), offset,
                                         direction = neighbor_cds.strand==origin_cds.strand,
                                         cds_name = neighbor_cds.cds_name,
                                         gene_name = neighbor_cds.gene_name)
                    gene2positions[neighbor_cds.gene_name].append(pos)
                row.append(pos)
            assert len(row) == 2 * self.DIST + 1
            matrix.append(row)

        self.origin_gene_name = origin_gene_name
        self.matrix = matrix
        self.shape = (len(self.matrix), 2 * self.DIST + 1)
        self.gene2positions = gene2positions

    def __repr__(self):
        return "<Matrix@{0}({1}x{2})>".format(self.origin_gene_name, self.shape[0], self.shape[1])

    def get_count_by_offset(self, offset):
        idx = offset + self.DIST
        return sum([(row[idx] is not None) for row in self.matrix])

    def get_positions_by_offset(self, offset):
        idx = offset + self.DIST
        return [row[idx] for row in self.matrix]

    def get_positions_by_gene_name(self, gene_name):
        return self.gene2positions[gene_name]

    def get_neighbor_gene_names(self):
        gene_name_set = set(self.gene2positions.keys())
        gene_name_set.remove(self.origin_gene_name)
        gene_name_set.remove(None) # default value for set_gene_name()
        return gene_name_set


def calc_bls(genome_names, tree):
    """
    calculate sum of branch length covered by a subset of leafs
    """

    genome_names = set(genome_names)
    if len(genome_names) <= 1:
        return 0

    subtree = tree.copy()
    subtree.prune(genome_names, preserve_branch_length=True)
    assert len(subtree) == len(genome_names)

    bls = 0
    stack = [subtree.get_tree_root()]
    while len(stack) > 0:
        parent = stack.pop()
        for child in parent.get_children():
            bls += subtree.get_distance(parent, child)
            stack.append(child)
    return bls

def set_split(cdss, split_fp):
    # TODO
    return cdss

def set_gene_name(cdss, ortho_fp):
    ortho_df = pd.read_csv(ortho_fp, sep='\t')
    cds2gene = dict([(cds, gene) for cds, gene in zip(ortho_df["cds_name"], ortho_df["gene_name"])])
    for cds in cdss:
        if cds.cds_name in cds2gene:
            cds.gene_name = cds2gene[cds.cds_name]
        else:
            cds.gene_name = None
    return cdss

def output_neighbor_df(edge_df, out_fp):
    with open(out_fp, 'w') as f:
        comment = ';'.join(["{}={}".format(key,val) for key,val in THRESH.items()])
        f.write("#{}\n".format(comment))
        edge_df.to_csv(f, index=False)
