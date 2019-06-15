CREATE TABLE hits (
        hit_id INTEGER NOT NULL,
        cds_id INTEGER,
        refseq_id INTEGER,
        length INTEGER,
        identity REAL,
        coverage REAL,
        PRIMARY KEY (hit_id)
);
