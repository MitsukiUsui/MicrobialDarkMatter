pragma journal_mode=memory;
pragma synchronous=0;
pragma cache_size=500000;

.read ./hits/schema.sql
.mode tabs
.import ./hits/data/hits.tsv hits
.read ./hits/index.sql

