import pandas as pd

MMSEQS_SCHEMA = [("qname", object), ("sname", object),
                 ("identity", float), ("length", int), ("mismatch", int), ("gap", int),
                 ("qstart", int), ("qend", int), ("sstart", int), ("send", int),
                 ("evalue", float), ("bitscore", int), ("qlength", int), ("slength", int)]


def read_mmseqs(fp, **kwargs):
    # trim schema based on the first line
    with open(fp, 'r') as f:
        line = f.readline().strip()
    schema = MMSEQS_SCHEMA[:len(line.split())]

    # overwrite default params
    params = {
        "header": None,
        "sep": '\t',
        "names": list(map(lambda x: x[0], schema)),
        "dtype": dict(schema)
    }
    for key, val in kwargs.items():
        params[key] = val

    df = pd.read_csv(fp, **params)
    return df
