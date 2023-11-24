make
./VAT makevatdb --dbtype prot --in ../data/protein/IPR034040_IPR034042.fa -d prot_db
/usr/bin/time -v ./VAT protein -d prot_db.vatf -q ../data/protein/IPR034042.fa -a match_prot_db -p 16 --id2 30
./VAT view -a match_prot_db.vatr -o match_prot
python3 ../scripts/protein_benchmark.py ../data/protein/IPR034042.fa match_prot
python3 ../scripts/overlap_blastp.py ../scripts/myresults.txt match_prot