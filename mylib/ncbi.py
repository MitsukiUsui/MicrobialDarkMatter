import pandas as pd
from logging import getLogger

from . import path

GENBANK_PATH = path.GENBANK_PATH
REFSEQ_PATH = path.REFSEQ_PATH
logger = getLogger(__name__)
genbankDAO = None
refseqDAO = None


class AssemblySummaryDAO:
    def __init__(self, assembly_summary_fp):
        self.df = pd.read_csv(assembly_summary_fp, sep='\t', skiprows=1, low_memory=False)
        self.acc2ftp = dict(self.df[["# assembly_accession", "ftp_path"]].values)
        logger.debug("loaded assembly summary from {}".format(assembly_summary_fp))

    def __contains__(self, accession):
        return accession in self.acc2ftp

    def build_ftp_filepath(self, accession, extension=None):
        possible_extension_set = {"fna", "faa"}

        # error handling
        if accession not in self.acc2ftp:
            logger.error("accession={} is not found".format(accession))
            return None
        if extension and extension not in possible_extension_set:
            logger.debug("extension={} is not allowed : {}".format(extension, possible_extension_set))
            return None

        ftp_direc = self.acc2ftp[accession]
        if extension is None:
            ftp_path = ftp_direc
        if extension == "fna":
            ftp_path = "{}/{}_genomic.fna.gz".format(ftp_direc, ftp_direc.split('/')[-1])
        elif extension == "faa":
            ftp_path = "{}/{}_protein.faa.gz".format(ftp_direc, ftp_direc.split('/')[-1])
        return ftp_path


def build_ftp_filepath(accession, extension=None):
    if genbankDAO is None or refseqDAO is None:
        __load()

    if accession in genbankDAO:
        return genbankDAO.build_ftp_filepath(accession, extension)
    elif accession in refseqDAO:
        return refseqDAO.build_ftp_filepath(accession, extension)
    else:
        logger.debug("accession={} not found".format(accession))
        return None


def __load():
    global genbankDAO
    global refseqDAO
    if genbankDAO is None:
        genbankDAO = AssemblySummaryDAO(GENBANK_PATH)
    if refseqDAO is None:
        refseqDAO = AssemblySummaryDAO(REFSEQ_PATH)
