import sys

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from . import path
sys.path.append(path.DB_DIREC)
from myschema import Project, Genome, Scaffold, Cds

DB_PATH=path.DB_PATH


def get_session():
    engine = create_engine('sqlite:///{}'.format(DB_PATH))
    Session = sessionmaker(bind=engine)
    return Session()

def get_connection():
    engine = create_engine('sqlite:///{}'.format(DB_PATH))
    return engine.connect()

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
