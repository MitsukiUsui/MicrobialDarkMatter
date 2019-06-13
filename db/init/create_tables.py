#!/usr/bin/env python3

import re
import sys
import pathlib
import logging

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

ROOT_PATH = pathlib.Path().joinpath('../../').resolve()
sys.path.append(str(ROOT_PATH))
import mylib
from mylib.db import Project, Genome, Scaffold, Cds
from mylib.db import IDManager
from mylib.path import build_local_filepath
from mylib.gff import parse_gff

LOGGER = logging.getLogger(__name__)
PID = IDManager("projects")
GID = IDManager("genomes")
SID = IDManager("scaffolds")
CID = IDManager("cdss")


def load_scaffolds(fna_fp):
    scaffolds = []
    for record in SeqIO.parse(fna_fp, 'fasta'):
        scaffolds.append(Scaffold(
            scaffold_id = SID.new(),
            genome_id = GID.get(),
            scaffold_name = record.id,
            length = len(record)
        ))
    return scaffolds

def load_scaffolds_faster(gff_fp):
    scaffolds = []
    pattern = re.compile(r"# Sequence Data: seqnum=[0-9]*;seqlen=([0-9]*);seqhdr=\"(.*)\"")
    with open(gff_fp, "r") as f:
        for line in f:
            m = re.match(pattern, line)
            if m is not None:
                scaffolds.append(Scaffold(
                    scaffold_id = SID.new(),
                    genome_id = GID.get(),
                    scaffold_name = m.group(2).split()[0],
                    length = int(m.group(1))
                ))
    return scaffolds

def load_cdss(gff_fp, scaffold2id):
    cdss = []
    for record in parse_gff(gff_fp):
        if record.type == "CDS":
            cdss.append(Cds(
                cds_id = CID.new(),
                genome_id = GID.get(),
                scaffold_id = scaffold2id[record.seqid],
                cds_name = record.attributes["cds_name"],
                start = record.start,
                end = record.end,
                length = record.end - record.start + 1,
                strand = record.strand
            ))
    return cdss

def append_records(records, out_fp):
    with open(out_fp, 'a') as f:
        for record in records:
            f.write('{}\n'.format(record))

def main(arg_fp, projects_fp, genomes_fp, scaffolds_fp, cdss_fp):
    arg_df = pd.read_csv(arg_fp, sep='\t', comment='#')
    LOGGER.info("found {} projects to load".format(len(arg_df)))
    for project_name, meta_fp in zip(arg_df["project_name"], arg_df["meta_fp"]):
        project = Project(project_id=PID.new(), project_name=project_name)
        meta_df = pd.read_csv(meta_fp, sep='\t')
        LOGGER.info("start {} with {} genomes".format(project_name, len(meta_df)))

        for _, genome_name in enumerate(tqdm(meta_df["genome_name"])):
            if _ >= 10:
                break
            LOGGER.debug("genome_name: {}".format(genome_name))
            genome = Genome(genome_id=GID.new(), project_id=PID.get(), genome_name=genome_name)
            fna_fp = build_local_filepath(genome_name, "fna").replace(".fna", ".dnaseq")
            gff_fp = build_local_filepath(genome_name, "gff")

            scaffolds = load_scaffolds(fna_fp)
            LOGGER.debug("loaded {} scaffolds from {}".format(len(scaffolds), fna_fp))
            scaffold2id = dict([(scaffold.scaffold_name, scaffold.scaffold_id) for scaffold in scaffolds])
            cdss = load_cdss(gff_fp, scaffold2id)
            LOGGER.debug("loaded {} cdss from {}".format(len(cdss), gff_fp))

            append_records(scaffolds, scaffolds_fp)
            append_records(cdss, cdss_fp)
            append_records([genome], genomes_fp)
        append_records([project], projects_fp)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, datefmt="%m/%d/%Y %I:%M:%S",
                            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    arg_fp = "./arg/create_tables.arg"
    projects_fp = "./data/projects.tsv"
    genomes_fp = "./data/genomes.tsv"
    scaffolds_fp = "./data/scaffolds.tsv"
    cdss_fp = "./data/cdss.tsv"
    main(arg_fp, projects_fp, genomes_fp, scaffolds_fp, cdss_fp)
