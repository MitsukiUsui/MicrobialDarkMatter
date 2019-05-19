import sys
from logging import getLogger

logger = getLogger(__name__)

DATA_DIREC = "/home/mitsuki/data/mag/genome"
DB_DIREC = "/home/mitsuki/mag/db"
DB_PATH = "{}/mag.db".format(DB_DIREC)

def build_local_filepath(genome_name, extension):
    possible_extension_set = set(["fna", "faa", "gff"])
    if extension in possible_extension_set:
        return "{0}/{1}/{1}.{2}".format(DATA_DIREC, genome_name, extension)
    else:
        logger.error("extension={} is not allowed : {}".format(extension, possible_extension_set))
        return None
