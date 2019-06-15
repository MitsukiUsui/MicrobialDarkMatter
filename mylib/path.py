import sys
from logging import getLogger

logger = getLogger(__name__)

DATA_DIREC = "/nfs_share/mitsuki/MicrobialDarkMatter/genome"
DB_LIB_DIREC = "/home/mitsuki/MicrobialDarkMatter/DB/init"
DB_PATH = "/home/mitsuki/MicrobialDarkMatter/DB/genome.db"
GENBANK_PATH = "/home/mitsuki/mag/fetch/data/raw/assembly_summary_genbank.txt"
REFSEQ_PATH = "/home/mitsuki/mag/fetch/data/raw/assembly_summary_refseq.txt"

def build_local_filepath(genome_name, extension=None):
    possible_extension_set = set(["fna", "faa", "gff"])

    local_direc = "{0}/{1}".format(DATA_DIREC, genome_name)
    if extension is None:
        return local_direc
    elif extension in possible_extension_set:
        return "{0}/{1}.{2}".format(local_direc, genome_name, extension)
    else:
        logger.error("extension={} is not allowed : {}".format(extension, possible_extension_set))
        return None
