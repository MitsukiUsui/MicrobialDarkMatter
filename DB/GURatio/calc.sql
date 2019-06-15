pragma journal_mode=memory;
pragma synchronous=0;
pragma cache_size=500000;

--------------------------------------------------------------------------------
-- cds_idとrefseq_idのマッチング (cdss.refseq_id)
--------------------------------------------------------------------------------
UPDATE cdss
SET refseq_id = (
    SELECT refseq_id
    FROM hits
    WHERE hits.cds_id = cdss.cds_id
    AND hits.coverage >= 0.2
);

--------------------------------------------------------------------------------
-- gkとfkの伝搬 (cdss.gk, cdss.fk)
--------------------------------------------------------------------------------
UPDATE cdss
SET
    gk = (
        SELECT gk
        FROM refseqs
        WHERE refseqs.refseq_id = cdss.refseq_id
    ),
    fk = (
        SELECT fk
        FROM refseqs
        WHERE refseqs.refseq_id = cdss.refseq_id
    )
WHERE refseq_id IS NOT NULL;

--------------------------------------------------------------------------------
-- scaffolds.cds_count, gk_count, fk_count
--------------------------------------------------------------------------------
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

--------------------------------------------------------------------------------
-- genomes.cds_count, gk_count, fk_count
--------------------------------------------------------------------------------
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

--------------------------------------------------------------------------------
-- clusters.cds_count, gk_count, fk_count
--------------------------------------------------------------------------------
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
