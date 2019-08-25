#!/usr/bin/env python3

import logging
import math
import pathlib
import random
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

ROOT_PATH = pathlib.Path().joinpath('../../').resolve()
sys.path.append(str(ROOT_PATH))
from mylib.db import load_genome_names_by_clade_name, load_cdss_by_genome_names
from splitlib import SegmentManager, Wcf

LOGGER = logging.getLogger(__name__)


def calc_loss(wcf1, wcf2):
    size = max(len(wcf1), len(wcf2))
    wcf1_ = np.ones(size)
    wcf2_ = np.ones(size)
    wcf1_[:len(wcf1)] = wcf1.to_array()
    wcf2_[:len(wcf2)] = wcf2.to_array()
    loss = abs(wcf1_ - wcf2_).sum()
    return loss


def fit(segment_manager, wcf_model):
    """
    update SegmentManager to follow wcf_model
    """

    total_member_count = segment_manager.get_member_count()
    used_member_count = 0
    for size in range(segment_manager.get_max_segment_size(), 5, -1):  # split segments in descending order
        segment_ids = segment_manager.get_segments_by_size(size)
        available_member_count = total_member_count * (1 - wcf_model[size - 1]) - used_member_count
        assert available_member_count >= 0
        keep_segment_count = math.floor(available_member_count / size)
        split_segment_count = len(segment_ids) - keep_segment_count

        used_member_count += size * keep_segment_count
        if split_segment_count > 0:
            for segment_id in random.sample(segment_ids, split_segment_count):
                idx = random.randint(3, size - 3)
                segment_manager.split(segment_id, idx)
    return segment_manager


def main(clade_name, model_fp, out_fp):
    genome_names = load_genome_names_by_clade_name(clade_name)
    LOGGER.info("loaded {} {} genomes".format(len(genome_names), clade_name))
    cdss = load_cdss_by_genome_names(genome_names)
    LOGGER.info("loaded {} cdss".format(len(cdss)))

    model_df = pd.read_csv(model_fp, sep='\t')
    wcf_model = Wcf(model_df["x"], model_df["y"])
    LOGGER.info("loaded model distribution from {}".format(model_fp))

    segment_manager = SegmentManager()
    scaffold2members = defaultdict(list)  # key: scaffold_id, val: list of cds_id
    for cds in cdss:
        scaffold2members[cds.scaffold_id].append(cds.cds_id)
    for _, members in scaffold2members.items():
        segment_manager.add(members)
    assert len(segment_manager) == len(scaffold2members)
    LOGGER.info("initialized segment manager with {} segments".format(len(segment_manager)))

    wcf_before = segment_manager.to_wcf()
    segment_manager = fit(segment_manager, wcf_model)
    wcf_after = segment_manager.to_wcf()
    LOGGER.info("updated sement manager to {} segments. Fitting loss: {} -> {}".format(
        len(segment_manager), calc_loss(wcf_model, wcf_before), calc_loss(wcf_model, wcf_after)))

    scaffold_df = pd.DataFrame(list(map(lambda cds: {"cds_id": cds.cds_id, "scaffold_id": cds.scaffold_id}, cdss)))
    segment_df = segment_manager.to_df().rename(columns={"member": "cds_id"})
    out_df = pd.merge(scaffold_df, segment_df)
    assert len(scaffold_df) == len(segment_df) == len(out_df)
    out_df.to_csv(out_fp, index=False, sep='\t')
    LOGGER.info("saved results to {}".format(out_fp))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, datefmt="%m/%d/%Y %I:%M:%S",
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    clade_name = "Enterobacterales"
    model_fp = "./splitdata/MGII.dist"
    out_fp = "./splitdata/00.map"
    main(clade_name, model_fp, out_fp)
