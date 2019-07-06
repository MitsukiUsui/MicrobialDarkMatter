import sys
from collections import defaultdict

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import pandas as pd

from . import path
sys.path.append(path.DB_LIB_DIREC)
from myschema import Project, Genome, Scaffold, Cds
DB_PATH = path.DB_PATH


def get_session(fp=None):
    engine = create_engine('sqlite:///{}'.format(fp if fp else DB_PATH))
    Session = sessionmaker(bind=engine)
    return Session()


def get_connection(fp=None):
    engine = create_engine('sqlite:///{}'.format(fp if fp else DB_PATH))
    return engine.connect()


class IDManager:
    """
    SQLITE3 ID Management utility class
    """

    def __init__(self, table_name):
        assert table_name in ("projects", "genomes", "scaffolds", "cdss", "hits", "refseqs")
        self.con = get_connection()
        self.table_name = table_name
        self.current_id = self._query_max_id()

    def __del__(self):
        self.con.close()

    def get(self):
        return self.current_id

    def new(self):
        self.current_id += 1
        return self.current_id

    def _query_max_id(self):
        query = "SELECT MAX({}_id) from {};".format(self.table_name[:-1], self.table_name)
        max_id = self.con.execute(query).fetchone()[0]
        max_id = max_id if max_id else 0  # for cases when table has no record
        return max_id


class CdsDAO:
    def __init__(self, cdss):
        self.cdss = cdss
        self.id2idx = defaultdict(lambda: None)
        self.name2idx = defaultdict(lambda: None)
        self.gene2idxs = defaultdict(list)
        for idx, cds in enumerate(self.cdss):
            self.id2idx[cds.cds_id] = idx
            self.name2idx[cds.cds_name] = idx
            if hasattr(cds, "gene_name"):
                self.gene2idxs[cds.gene_name].append(idx)

    def get_cds_by_idx(self, idx):
        if isinstance(idx, int) and 0 <= idx < len(self.cdss):
            return self.cdss[idx]
        else:
            return None

    def get_cds_by_cds_id(self, cds_id):
        return self.get_cds_by_idx(self.id2idx[cds_id])

    def get_cds_by_cds_name(self, cds_name):
        return self.get_cds_by_idx(self.name2idx[cds_name])

    def get_cdss_by_gene_name(self, gene_name):
        return list(map(lambda idx: self.get_cds_by_idx(idx), self.gene2idxs[gene_name]))

    def get_neighbor_cds(self, origin_cds, offset):
        neighbor_cds_id = origin_cds.cds_id + offset if origin_cds.strand == '+' else origin_cds.cds_id - offset
        neighbor_cds = self.get_cds_by_cds_id(neighbor_cds_id)
        if neighbor_cds is not None and neighbor_cds.scaffold_id == origin_cds.scaffold_id:
            return neighbor_cds
        else:
            return None


def load_name2id(table_name, default=-1, con=None):
    assert table_name in ("projects", "genomes", "scaffolds", "cdss", "refseqs")
    create_tmp_con = con is None
    if create_tmp_con:
        con = get_connection()

    col_id = "{}_id".format(table_name[:-1])
    col_name = "{}_name".format(table_name[:-1])
    query = "SELECT {}, {} FROM {};".format(col_id, col_name, table_name)
    df = pd.read_sql_query(query, con)
    name2id = defaultdict(lambda: default)
    for name, id_ in zip(df[col_name], df[col_id]):
        name2id[name] = id_

    if create_tmp_con:
        con.close()
    return name2id


def load_genome_names_by_clade_name(clade_name, con=None):
    create_tmp_con = con is None
    if create_tmp_con:
        con = get_connection()

    query = "SELECT genome_name FROM clades WHERE clade_name='{}'".format(clade_name)
    clades_df = pd.read_sql_query(query, con)
    genome_names = list(clades_df["genome_name"])

    if create_tmp_con:
        con.close()
    return genome_names


def load_genomes_by_genome_names(genome_names, session=None):
    create_tmp_session = session is None
    if create_tmp_session:
        session = get_session()

    genomes = session.query(Genome).filter(Genome.genome_name.in_(genome_names)).all()
    assert len(genomes) == len(genome_names)

    if create_tmp_session:
        session.close()
    return genomes


def load_cdss_by_genome_names(genome_names, session=None):
    create_tmp_session = session is None
    if create_tmp_session:
        session = get_session()

    genomes = load_genomes_by_genome_names(genome_names, session)
    genome_ids = [genome.genome_id for genome in genomes]
    cdss = session.query(Cds).filter(Cds.genome_id.in_(genome_ids)).all()

    if create_tmp_session:
        session.close()
    return cdss
