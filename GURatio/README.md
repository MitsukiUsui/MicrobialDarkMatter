cat ./arg/refseqs.arg |xargs cat > ./data/refseqs.faa
seqkit fx2tab ./data/refseqs.faa -nl > ./data/refseqs.stat

