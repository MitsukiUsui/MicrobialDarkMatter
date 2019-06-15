from sqlalchemy import create_engine
from sqlalchemy import Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()


class Project(Base):
    __tablename__ = "projects"
    project_id = Column(Integer, primary_key=True)
    project_name = Column(String)

    def __init__(self, project_id, project_name):
        self.project_id = project_id
        self.project_name = project_name

    def __str__(self):
        return '\t'.join(map(str, [self.project_id, self.project_name]))

    def __repr__(self):
        return "<Project(name={})>".format(self.project_name)

class Genome(Base):
    __tablename__ = "genomes"
    genome_id = Column(Integer, primary_key=True)
    project_id = Column(Integer)
    genome_name = Column(String)

    def __init__(self, genome_id, project_id, genome_name):
        self.genome_id = genome_id
        self.project_id = project_id
        self.genome_name = genome_name

    def __str__(self):
        return '\t'.join(map(str, [self.genome_id, self.project_id, self.genome_name]))

    def __repr__(self):
        return "<Genome(name={})>".format(self.genome_name)

class Scaffold(Base):
    __tablename__ = "scaffolds"
    scaffold_id = Column(Integer, primary_key=True)
    genome_id = Column(Integer)
    scaffold_name = Column(String)
    length = Column(Integer)

    def __init__(self, scaffold_id, genome_id, scaffold_name, length):
        self.scaffold_id = scaffold_id
        self.genome_id = genome_id
        self.scaffold_name = scaffold_name
        self.length = length

    def __str__(self):
        return '\t'.join(map(str, [self.scaffold_id, self.genome_id, self.scaffold_name, self.length]))

    def __repr__(self):
        return "<Scaffold(name={})>".format(self.scaffold_name)

class Cds(Base):
    __tablename__ = "cdss"
    cds_id = Column(Integer, primary_key=True)
    genome_id = Column(Integer)
    scaffold_id = Column(Integer)
    cds_name = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    length = Column(Integer)
    strand = Column(String)

    def __init__(self, cds_id, genome_id, scaffold_id, cds_name, start, end, length, strand):
        self.cds_id = cds_id
        self.genome_id = genome_id
        self.scaffold_id = scaffold_id
        self.cds_name = cds_name
        self.start = start
        self.end = end
        self.length = length
        self.strand = strand

    def __str__(self):
        return '\t'.join(map(str, [self.cds_id, self.genome_id, self.scaffold_id, self.cds_name,
                                    self.start, self.end, self.length, self.strand]))

    def __repr__(self):
        return "<Cds(name={})>".format(self.cds_name)
