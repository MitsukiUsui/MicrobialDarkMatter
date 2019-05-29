CREATE TABLE projects (
    project_id INTEGER NOT NULL,
    project_name TEXT NOT NULL UNIQUE,
    PRIMARY KEY (project_id)
);

CREATE TABLE genomes(
    genome_id INTEGER NOT NULL,
    project_id INTEGER NOT NULL,
    genome_name TEXT NOT NULL UNIQUE,
    PRIMARY KEY (genome_id)
);

CREATE TABLE scaffolds(
    scaffold_id INTEGER NOT NULL,
    genome_id INTEGER NOT NULL,
    scaffold_name TEXT NOT NULL UNIQUE,
    length INTEGER,
    PRIMARY KEY (scaffold_id)
);

CREATE TABLE cdss(
    cds_id INTEGER NOT NULL,
    genome_id INTEGER NOT NULL,
    scaffold_id INTEGER NOT NULL,
    cds_name TEXT NOT NULL UNIQUE,
    start INTEGER,
    end INTEGER,
    length INTEGER,
    strand TEXT,
    PRIMARY KEY (cds_id)
);
