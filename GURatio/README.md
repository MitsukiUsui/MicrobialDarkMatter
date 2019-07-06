Create reference .faa from RefSeq proteins
```
./list_refseqs.py ./arg/refseqs.arg
cat ./arg/refseqs.arg |xargs cat > ./data/refseqs.faa
seqkit fx2tab ./data/refseqs.faa -nl > ./data/refseqs.stat
```

Create query .faa
```
./list_queries.py ./arg/queries.arg
cat ./arg/queries.arg |xargs cat > ./data/queries.faa
```