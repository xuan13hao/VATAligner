# Introduction

## DNA
./VAT makevatdb --dbtype nucl --in ../data/test_all.fa -d mt
./VAT dna -d mt.vatf -q ../data/test_forward_reads.fa -a match
./VAT view -a match.vatr -o match
vim match

## DNA Zip files
./VAT makevatdb --dbtype nucl --in ../data/ref.fa.gz -d mt
./VAT dna -d mt.vatf -q ../data/query.fa.gz -a match -p 4
./VAT view -a match.vatr -o match
vim match

## Whole-genome alignment
./VAT makevatdb --dbtype nucl --in ../data/test_all.fa -d mt
./VAT dna -d mt.vatf -q ../data/test_forward_reads.fa -a match --whole-genome
./VAT view -a match.vatr -o match
vim match

## Protein
./VAT makevatdb --dbtype prot --in ../data/pfam_ref.fa -d pfam
./VAT protein -d pfam.vatf -q ../data/pfam_test.fa -a match -p 4
./VAT view -a match.vatr -o match
vim match

