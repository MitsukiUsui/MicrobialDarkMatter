import pandas as pd
from logging import getLogger

from . import path

GENBANK_PATH=path.GENBANK_PATH
logger = getLogger(__name__)
genbankDAO = None


class GenbankDAO:
    def __init__(self):
        self.df = pd.read_csv(GENBANK_PATH, sep='\t', skiprows=1)
        self.acc2ftp = dict(self.df[["# assembly_accession", "ftp_path"]].values)
        logger.debug("loaded genbank data from {}".format(GENBANK_PATH))

    def build_ftp_filepath(self, accession, extension=None):
        possible_extension_set = set(["fna", "faa"])

        # error handling
        if accession not in self.acc2ftp:
            logger.error("accession={} is not found".format(accession))
            return None
        if extension and extension not in possible_extension_set:
            logger.error("extension={} is not allowed : {}".format(extension, possible_extension_set))
            return None

        ftp_direc = self.acc2ftp[accession]
        if extension is None:
            ftp_path = ftp_direc
        if extension == "fna":
            ftp_path  = "{}/{}_genomic.fna.gz".format(ftp_direc, ftp_direc.split('/')[-1])
        elif extension == "faa":
            ftp_path  = "{}/{}_protein.faa.gz".format(ftp_direc, ftp_direc.split('/')[-1])
        return ftp_path

def build_ftp_filepath(accession, extension=None):
    if genbankDAO is None:
        __load()
    return genbankDAO.build_ftp_filepath(accession, extension)

def __load():
    global genbankDAO
    if genbankDAO is None:
        genbankDAO = GenbankDAO()
