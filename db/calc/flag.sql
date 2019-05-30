pragma journal_mode=memory;
pragma synchronous=0;
pragma cache_size=500000;

--------------------------------------------------------------------------------
-- cds_idとrefseq_idのマッチング (cdss.refseq_id)
--------------------------------------------------------------------------------
ALTER TABLE cdss ADD COLUMN refseq_id INTEGER;
UPDATE cdss SET refseq_id = NULL;

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
ALTER TABLE cdss ADD COLUMN gk INTEGER DEFAULT 0;
ALTER TABLE cdss ADD COLUMN fk INTEGER DEFAULT 0;
UPDATE cdss SET gk = 0;
UPDATE cdss SET fk = 0;

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

