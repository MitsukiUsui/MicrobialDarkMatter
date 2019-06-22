import sys
from collections import defaultdict

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import pandas as pd

from . import path
sys.path.append(path.DB_LIB_DIREC)
from myschema import Project, Genome, Scaffold, Cds

DB_PATH=path.DB_PATH


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
        max_id = max_id if max_id else 0 # for cases when table has no record
        return max_id

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

def load_genomes(genome_names, session=None, to_dict=False):
    create_tmp_session = session is None
    if create_tmp_session:
        session = get_session()

    genomes = session.query(Genome).filter(Genome.genome_name.in_(genome_names)).all()
    assert len(genomes) == len(genome_names)
    if to_dict:
        return dict([(genome.genome_id, genome) for genome in genomes])

    if create_tmp_session:
        session.close()
    else:
        return genomes

def load_cdss_by_genome_names(genome_names, session=None, to_dict=False):
    create_tmp_session = session is None
    if create_tmp_session:
        session = get_session()

    genomes = load_genomes(genome_names, session)
    genome_ids = [genome.genome_id for genome in genomes]
    cdss = session.query(Cds).filter(Cds.genome_id.in_(genome_ids)).all()
    if to_dict:
        return dict([(cds.cds_id, cds) for cds in cdss])

    if create_tmp_session:
        session.close()
    else:
        return cdss
