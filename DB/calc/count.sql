pragma journal_mode=memory;
pragma synchronous=0;
pragma cache_size=500000;

--scaffolds.cds_count, gk_count, fk_count
ALTER TABLE scaffolds ADD COLUMN cds_count INTEGER;
ALTER TABLE scaffolds ADD COLUMN gk_count INTEGER;
ALTER TABLE scaffolds ADD COLUMN fk_count INTEGER;

UPDATE scaffolds
SET
    cds_count = (
        SELECT COUNT(cds_id)
        FROM cdss
        WHERE cdss.scaffold_id = scaffolds.scaffold_id
    ),
    gk_count = (
        SELECT SUM(gk)
        FROM cdss
        WHERE cdss.scaffold_id = scaffolds.scaffold_id
    ),
    fk_count = (
        SELECT SUM(fk)
        FROM cdss
        WHERE cdss.scaffold_id = scaffolds.scaffold_id
    )
;


--genomes.cds_count, gk_count, fk_count
ALTER TABLE genomes ADD COLUMN cds_count INTEGER;
ALTER TABLE genomes ADD COLUMN gk_count INTEGER;
ALTER TABLE genomes ADD COLUMN fk_count INTEGER;

UPDATE genomes
SET
    cds_count = (
        SELECT COUNT(cds_id)
        FROM cdss
        WHERE cdss.genome_id = genomes.genome_id
    ),
    gk_count = (
        SELECT SUM(gk)
        FROM cdss
        WHERE cdss.genome_id = genomes.genome_id
    ),
    fk_count = (
        SELECT SUM(fk)
        FROM cdss
        WHERE cdss.genome_id = genomes.genome_id
    )
;

--clusters.cds_count, gk_count, fk_count
ALTER TABLE clusters ADD COLUMN cds_count INTEGER;
ALTER TABLE clusters ADD COLUMN gk_count INTEGER;
ALTER TABLE clusters ADD COLUMN fk_count INTEGER;

UPDATE clusters
SET
    cds_count = (
        SELECT COUNT(cds_id)
        FROM cdss
        WHERE cdss.cluster_id = clusters.cluster_id
    ),
    gk_count = (
        SELECT SUM(gk)
        FROM cdss
        WHERE cdss.cluster_id = clusters.cluster_id
    ),
    fk_count = (
        SELECT SUM(fk)
        FROM cdss
        WHERE cdss.cluster_id = clusters.cluster_id
    )
;
