import pandas as pd

def to_hits_df(fp):
    df = pd.read_csv(fp, header=None, sep='\t',
                     names=["qname", "sname", "identity", "length", "mismatch", "gap",
                            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlength", "slength"])
    return df
