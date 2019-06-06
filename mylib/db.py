import sys

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

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
        assert table_name in ("projects", "genomes", "scaffolds", "cdss", "hits")
        self.current_id = query_max_id(table_name)
    def get(self):
        return self.current_id
    def new(self):
        self.current_id += 1
        return self.current_id

def query_max_id(table_name, con=None):
    create_tmp_con = con is None

    if create_tmp_con:
        con = get_connection()
    query = "SELECT MAX({}_id) from {};".format(table_name[:-1], table_name)
    max_id = con.execute(query).fetchone()[0]
    max_id = max_id if max_id else 0 # for cases when table has no record
    if create_tmp_con:
        con.close()

    return max_id

def load_genomes(genome_names, session=None, to_dict=False):
    create_tmp_session = session is None

    if create_tmp_session:
        session = get_session()
    genomes = session.query(Genome).filter(Genome.genome_name.in_(genome_names)).all()
    assert len(genomes) == len(genome_names)
    if create_tmp_session:
        session.close()

    if to_dict:
        return dict([(genome.genome_id, genome) for genome in genomes])
    else:
        return genomes

def load_cdss_by_genome_names(genome_names, session=None, to_dict=False):
    create_tmp_session = session is None

    if create_tmp_session:
        session = get_session()
    genomes = load_genomes(genome_names, session)
    genome_ids = [genome.genome_id for genome in genomes]
    cdss = session.query(Cds).filter(Cds.genome_id.in_(genome_ids)).all()
    if create_tmp_session:
        session.close()

    if to_dict:
        return dict([(cds.cds_id, cds) for cds in cdss])
    else:
        return cdss
