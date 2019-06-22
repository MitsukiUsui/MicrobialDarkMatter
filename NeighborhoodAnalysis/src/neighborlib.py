#!/usr/bin/env python3

from logging import getLogger
from collections import defaultdict, Counter

import numpy as np

THRESH = {
    "SCORE" : 0.8,
    "SIZE" : 10,
    "DIST" : 5
}
LOGGER = getLogger(__name__)


class MatrixPosition:
    """
    store position in NeighborhoodMatrix (NxM)
    """

    def __init__(self, i, j, offset, direction):
        self.i = i
        self.j = j
        self.offset = offset
        self.direction = direction

    def __lt__(self, other):
        assert self.i == other.i
        return abs(self.offset) < abs(other.offset)

    def __repr__(self):
        return "<MatrixPosition({},{})>".format(self.i, self.j)

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

def is_neighborhood(cds1, cds2):
    if cds1.scaffold_id != cds2.scaffold_id:
        return False
    if abs(cds1.cds_id - cds2.cds_id) > THRESH["DIST"]:
        return False
    return True

def scan_neighborhoods(origin_gene_name, cdsDAO):
    """
    scan neighborhoods and create observation matrix and inverted indexes
    """

    origin_cdss = cdsDAO.get_cdss_by_gene_name(origin_gene_name)
    H = len(origin_cdss)
    W = THRESH["DIST"] * 2 + 1
    obs_mat = np.zeros((H, W)).astype(bool) #initialize by False
    gene2positions = defaultdict(list) #key: gene_name, val: list of neighborhood matrix positions

    for i in range(H):
        for j in range(W):
            origin_cds = origin_cdss[i]
            offset = j - THRESH["DIST"]
            neighbor_cds_id = origin_cds_id + offset if origin_cds.strand == '+' else origin_cds_id - offset
            neighbor_cds = cdsDAO.get_cds_by_cds_id(neighbor_cds_id)
            if not(neighbor_cds is None) and is_neighborhood(origin_cds, neighbor_cds):
                obs_mat[i][j] = True
                if hasattr(neighbor_cds, "gene_name"):
                    pos = MatrixPosition(i, j, offset, origin_cds.strand==neighbor_cds.strand)
                    gene2positions[neighbor_cds.gene_name].append(pos)
    return obs_mat, gene2positions

def find_most_common_position(positions):
    def get_relationship(offset, direction):
        if direction == True:
            return "cooriented"
        elif offset < 0:
            return "divergent"
        else:
            return "convergent"

    assert len(positions) > 0
    c = Counter([(pos.offset, pos.direction) for pos in positions])
    (offset, direction), count = c.most_common()[0]

    return {
        "top_offset": offset,
        "top_relationship": get_relationship(offset, direction),
        "top_ratio": count / len(positions)
    }

def set_split(cdss, split_fp):
    # TODO
    return cdss

def set_gene_name(cdss, ortho_fp):
    # TODO
    return cdss

def output_neighbor_df(edge_df, out_fp):
    with open(out_fp, 'w') as f:
        comment = ';'.join(["{}={}".format(key,val) for key,val in THRESH.items()])
        f.write("#{}\n".format(comment))
        edge_df.to_csv(f, index=False)
