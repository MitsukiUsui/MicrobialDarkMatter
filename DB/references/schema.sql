CREATE TABLE refseqs (
    refseq_id INTEGER NOT NULL,
    refseq_name TEXT NOT NULL UNIQUE,
    description TEXT,
    lca TEXT,
    gk INTEGER,
    fk INTEGER,
    PRIMARY KEY (refseq_id)
);
