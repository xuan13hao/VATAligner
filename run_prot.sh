./VAT makevatdb --dbtype prot --in ../data/protein/IPR034040_IPR034042.fa -d prot_db
/usr/bin/time -v ./VAT protein -d prot_db.vatf -q ../data/protein/IPR034042.fa -a match_prot_db -p 16 --top 98 -t /dev/shm --max-target-seqs 4000
./VAT view -a match_prot_db.vatr -o match_prot
