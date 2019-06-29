#!/usr/bin/env python3

import sys
import pandas as pd
import networkx as nx

class Edge:
    def __init__(self, s, t):
        self.s = s
        self.t = t
    def __eq__(self, other):
        return self.s == other.s and self.t == other.t
    def __repr__(self):
        return "<Edge({}->{})>".format(self.s, self.t)
    def __hash__(self):
        return hash(self.s+self.t) #assume s+t != t+s ?
    def __contains__(self, item):
        return self.s==item or self.t==item
    def rev(self):
        return Edge(self.t, self.s)


def filter_neighbor_df(neighbor_df, score, bls, only_bi):
    msk = (neighbor_df["score"] >= score) & (neighbor_df["bls"] >= bls)
    edges = set([Edge(row["x"], row["y"]) for _, row in neighbor_df[msk].iterrows()])
    if only_bi:
        biedges = set([e for e in edges if e.rev() in edges])
        msk = [Edge(row["x"], row["y"]) in biedges for _, row in neighbor_df.iterrows()]
    else:
        msk = [Edge(row["x"], row["y"]) in edges for _, row in neighbor_df.iterrows()]
    return neighbor_df[msk]


def main(neighbor_fp, meta_fp, edge_fp, score, bls, only_bi):
    print("params: score={}, bls={}, bidirec={}".format(score, bls, only_bi))
    neighbor_df = pd.read_csv(neighbor_fp, comment="#")
    filtered_df = filter_neighbor_df(neighbor_df, score=score, bls=bls, only_bi=only_bi)
    print("filtered to {}%".format(len(filtered_df)/len(neighbor_df)*100))

    meta_df = pd.read_csv(meta_fp, sep='\t')
    meta_df = meta_df[["x", "y", "is_match", "x_path_desc","y_path_desc", "common_path_desc"]].fillna('')
    join_df = pd.merge(filtered_df, meta_df, on=["x", "y"])
    assert len(join_df) == len(filtered_df)

    if only_bi:
        G = nx.Graph()
    else:
        G = nx.DiGraph()
    for _, row in join_df.iterrows():
        G.add_edge(row["x"], row["y"],is_match=row["is_match"],
                   x_path_desc=row["x_path_desc"], y_path_desc=row["y_path_desc"], common_path_desc=row["common_path_desc"],
                   bls=row["bls"], score=row["score"], total=row["total"], found=row["found"])

    nx.write_graphml(G, edge_fp)
    print("DONE: output {}".format(edge_fp))


if __name__=="__main__":
    neighbor_fp = sys.argv[1]
    meta_fp = sys.argv[2]
    edge_fp = sys.argv[3]
    main(neighbor_fp, meta_fp, edge_fp, score=0.3, bls=-1, only_bi=False)

