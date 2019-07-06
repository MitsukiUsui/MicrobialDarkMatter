#!/usr/bin/env python3

from logging import getLogger
from collections import defaultdict

import numpy as np
import pandas as pd

LOGGER = getLogger(__name__)


class MatrixPosition:
    def __init__(self, i, j, offset, direction, origin_name, cds_name, gene_name):
        self.i = i
        self.j = j
        self.offset = offset
        self.direction = direction  # ToDo rename to is_forward
        self.origin_name = origin_name
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

        self.matrix: 2D (HxW) matrix of MatrixPosition (or None), where H = number of origin cdss and W = 2 * DIST + 1
        self.indicator_matrix: numpy 2D (HxW) matrix , where -1 for missing, 0 for observed. This is the template for to_indicator_matrix(self, gene_name)
        self.offset2count: key: offset, val: observed count
        self.offset2positions: key: offset, val: list of MatrixPosition and None
        self.gene2positions: key: gene_name, val: list of MatrixPosition
        """

        self.origin_gene_name = origin_gene_name

        # initialize matrix
        matrix = []
        for i, origin_cds in enumerate(cdsDAO.get_cdss_by_gene_name(origin_gene_name)):
            row = []
            for j, offset in enumerate(range(-self.DIST, self.DIST+1)):
                neighbor_cds = cdsDAO.get_neighbor_cds(origin_cds, offset)
                if neighbor_cds is None:
                    pos = None
                else:
                    pos = MatrixPosition(i, j, offset,
                                         direction = neighbor_cds.strand==origin_cds.strand,
                                         origin_name = origin_cds.cds_name,
                                         cds_name = neighbor_cds.cds_name,
                                         gene_name = neighbor_cds.gene_name)
                row.append(pos)
            assert len(row) == 2 * self.DIST + 1
            matrix.append(row)
        self.matrix = matrix
        self.shape = (len(self.matrix), 2 * self.DIST + 1)

        # initialize index structures
        self.indicator_matrix = -np.ones((self.shape[0], self.shape[1])).astype(float)
        self.offset2count = defaultdict(lambda:0)
        self.offset2positions = defaultdict(list)
        self.gene2positions = defaultdict(list)
        for i, row in enumerate(self.matrix):
            for j, pos in enumerate(row):
                offset = j - self.DIST
                self.offset2positions[offset].append(pos)
                if pos is not None:
                    self.indicator_matrix[i][j] = 0
                    self.offset2count[offset] += 1
                    self.gene2positions[pos.gene_name].append(pos)

    def __repr__(self):
        return "<Matrix@{0}({1}x{2})>".format(self.origin_gene_name, self.shape[0], self.shape[1])

    def get_count_by_offset(self, offset):
        """
        :param offset: should between [-DIST, DIST]
        :return: number of cdss found at the offset
        """

        return self.offset2count[offset]

    def get_positions_by_offset(self, offset, dropna=False):
        """
        :param offset:
        :param dropna: remove None if true
        :return: all MatrixPositions or None at the offset
        """

        if dropna:
            return list(filter(lambda pos: pos is not None, self.offset2positions[offset]))
        else:
            return self.offset2positions[offset]

    def get_positions_by_gene_name(self, gene_name):
        """
        :param gene_name:
        :return: all MatrixPoisitions where the gene found
        """

        return self.gene2positions[gene_name]

    def get_neighbor_gene_names(self):
        """
        list all the neighborhood gene names which appeared at least once in the matrix.
        """

        gene_name_set = set(self.gene2positions.keys())
        gene_name_set.remove(self.origin_gene_name)
        if None in gene_name_set:  # default value for set_gene_name()
            gene_name_set.remove(None)
        return gene_name_set

    def to_indicator_matrix(self, gene_name):
        """
        convert to indicator matrix, where 1 stands for target gene, 0 for other gene, and -1 for missing
        """

        ret = self.indicator_matrix.copy()
        for pos in self.get_positions_by_gene_name(gene_name):
            ret[pos.i][pos.j] = 1
        return ret

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
