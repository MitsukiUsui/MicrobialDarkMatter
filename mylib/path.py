from logging import getLogger

LOGGER = getLogger(__name__)
GENOME_DIREC = "/nfs_share/mitsuki/MicrobialDarkMatter/genome"
CLADE_DIREC = "/nfs_share/mitsuki/MicrobialDarkMatter/clade"
DB_LIB_DIREC = "/home/mitsuki/MicrobialDarkMatter/DB/init"
DB_PATH = "/home/mitsuki/MicrobialDarkMatter/DB/genome.db"
GENBANK_PATH = "/home/mitsuki/mag/fetch/data/raw/assembly_summary_genbank.txt"
REFSEQ_PATH = "/home/mitsuki/mag/fetch/data/raw/assembly_summary_refseq.txt"


def build_local_filepath(genome_name, extension=None):
    possible_extension_set = {"fna", "faa", "gff"}

    local_direc = "{0}/{1}".format(GENOME_DIREC, genome_name)
    if extension is None:
        return local_direc
    else:
        assert extension in possible_extension_set
        return "{0}/{1}.{2}".format(local_direc, genome_name, extension)


def build_clade_filepath(clade_name):
    local_direc = "{0}/{1}".format(CLADE_DIREC, clade_name)
    return local_direc
