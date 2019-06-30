#!/usr/bin/env python3

import sys
import pathlib
import logging
import argparse
from collections import Counter

import pandas as pd

ROOT_PATH = pathlib.Path().joinpath('../../').resolve()
sys.path.append(str(ROOT_PATH))
from mylib.db import load_genome_names_by_clade_name, load_cdss_by_genome_names

LOGGER = logging.getLogger(__name__)

class SegmentManager:
    """
    * それぞれの要素がどのsegmentに属するか記憶する（self.segment2members）
    * サイズによってsegmentがひける（self.size2segments）
    * segmentの追加、削除、さらに分割が可能
    """

    def __init__(self):
        self.next_id = 0
        self.segment2members = dict() #key: segment_id, val: list of members
        self.size2segments = defaultdict(set) #key: size, val: set of segment_ids of the size
        self.total_member_count = 0

	def __len__(self):
		return len(self.segment2members)

    def add(self, members):
        segment_id = self.next_id
		size = len(members)
        self.next_id += 1
        self.segment2members[segment_id] = members
        self.size2segments[size].add(segment_id)
        self.total_member_count += size
        return segment_id

    def delete(self, segment_id):
        size = len(self.segment2members[segment_id])
        del self.segment2members[segment_id]
        self.size2segments[size].remove(segment_id)
        self.total_member_count -= size

    def split(self, segment_id, idx):
        members = self.get_members_by_id(segment_id)
        self.delete(segment_id)
        self.add(members[:idx])
        self.add(members[idx:])

    def get_members_by_id(self, segment_id):
        return self.segment2members[segment_id]

	def get_segments_by_size(self, size):
        return self.size2segments[size]

    def get_max_segment_size(self):
        return max(self.size2segments.keys())

    def get_total_member_count(self):
        return self.total_member_count

    def get_wcf(self):
        x = []
        y = []
        for size, segment_ids in self.size2segments.items():
            if len(segment_ids) > 0:
                x.append(size)
                y.append(len(segment_ids))
		wcf = calc_wcf(x, y)
        return wcf

    def to_df(self):
        if len(self.segment2members) == 0:
            return pd.DataFrame(names=["segment_id", "member"])

        records = []
        for segment_id, members in self.segment2members.items():
            records += [{"segment_id": segment_id, "member": m} for m in members]
        df = pd.DataFrame(records, names=["segment_id", "member"])
        return df

def calc_wcf(x, y):
    """
    calculate weighted cumlative frequencies
    """

    assert len(x) == len(y)
    weights = [x*y for x, y in zip(x, y)]
    res = stats.cumfreq(x, numbins=max(x)+1, defaultreallimits=(0, max(x)+1), weights=weights)
    wcf = res.cumcount / res.cumcount[-1]
    return wcf

def calc_loss(wcf1, wcf2):
    size = max(len(wcf1), len(wcf2))
    wcf1_ = np.ones(size)
    wcf2_ = np.ones(size)
    wcf1_[:len(wcf1)] = wcf1
    wcf2_[:len(wcf2)] = wcf2
    loss = abs(wcf1_ - wcf2_).sum()
    return loss

def fit(segment_manager, wcf_model):
	"""
	update sement manager to follow model_wcf
	"""

	total_member_count = segment_manager.get_total_member_count()
	used_member_count = 0
	for size in reversed(6, segment_manager.get_max_segment_size()): # split segments in descending order
		segment_ids = segment_manager.get_segments_by_size(size)
		if size > len(wcf_model):
			available_member_count = 0
		else:
			available_member_count = max(0, total_member_count * (1 - wcf_model[size-1]) - used_member_count)
		keep_segment_count = math.floor(available_member_count / size)
		split_segment_count = len(segment_ids) - keep_segment_count

		used_member_count += size * keep_segment_count
		if split_segment_count > 0:
            for segment_id in random.sample(segment_ids, split_segment_count):
				idx = random.randint(3, size-3)
                segment_manager.split(segment_id, idx)
	return segment_manager

def main(clade_name, model_fp, split_fp):
	genome_names = load_genome_names_by_clade_name(clade_name)
    LOGGER.info("loaded {} {} genomes".format(len(genome_names), clade_name))
    cdss = load_cdss_by_genome_names(genome_names)
    LOGGER.info("loaded {} cdss".format(len(cdss)))

    #load model distribution to fit
    model_df = pd.read_csv(model_fp, names=["x","y"])
    wcf_model = calc_wcf(model_df["x"], model_df["y"])
    LOGGER.info("loaded model distribution from {}".format(model_fp))

    #initialize SegmentManager with scaffolds as seed
    segment_manager = SegmentManager()
    scaffold2members = defaultdict(list) #key: scaffold_id, val: list of cds_id
    for cds in cdss:
        scaffold2members[cds.scaffold_id].append(cds.cds_id)
    for _, members in scaffold2members.items():
        segment_manager.add(members)
	assert len(segment_manager) == len(scaffold2members)
    LOGGER.info("initialized segment manager with {} segments".format(len(segment_manager)))

    #split segments
    wcf_before = segment_manager.get_wcf()
    segment_manager = fit(segment_manager, wcf_model)
    wcf_after = segment_manager.get_wcf()
    LOGGER.info("updated sement manager to {} segments. Fitting loss = {}".format(len(segment_manager), calc_loss(wcf_model, wcf_after))

    #output map to split_fp
    scaffold_df = pd.DataFrame(map(lambda cds: {"cds_id":cds.cds_id, "scaffold_id":cds.scaffold_id}), cdss)
    segment_df = segment_manager.to_df().rename(columns={"member":"cds_id"})
    out_df = pd.merge(scaffold_df, segment_df)
    assert len(scaffold_df) == len(segment_df) == len(out_df)
    out_df.to_csv(split_fp, index=False, sep='\t')
    LOGGER.info("saved results to {}".format(split_fp))

if __name__ == "__main__":
	logging.basicConfig(level=logging.INFO, datefmt="%m/%d/%Y %I:%M:%S",
                            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
	clade_name = ""
    model_fp = "./split/MGII_95.dist"
    out_fp = "./split/MGII_95.csv"
    main(clade_name, model_fp, split_fp)

