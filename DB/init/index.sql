CREATE INDEX ix_genomes_genome_name ON genomes (genome_name);
CREATE INDEX ix_scaffolds_scaffold_name ON scaffolds (scaffold_name);
CREATE INDEX ix_scaffolds_genome_id ON scaffolds (genome_id);
CREATE INDEX ix_cdss_cds_name ON cdss (cds_name);
CREATE INDEX ix_cdss_genome_id ON cdss (genome_id);
CREATE INDEX ix_cdss_scaffold_id ON cdss (scaffold_id);
