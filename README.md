# VAT Aligners

We propose to develop a novel multi-purpose aligner that is easy to use and has matched efficiency and alignment performance w.r.t its single-purpose competitors.


## Install
    1. cd src/
    2. install boost: ./pre-build_boost

## DNA
```console
VAT makevatdb --dbtype nucl --in ../data/test_all.fa -d mt
VAT dna -d mt.vatf -q ../data/test_forward_reads.fa -a match
VAT view -a match.vatr -o match
vim match
```
## DNA for multiple threads
```console
VAT makevatdb --dbtype nucl --in ../data/ref.fa.gz -d mt
VAT dna -d mt.vatf -q ../data/query.fa.gz -a match -p 4
VAT view -a match.vatr -o match
vim match
```
## Whole-genome alignment
```console
VAT makevatdb --dbtype nucl --in ../data/test_all.fa -d mt
VAT dna -d mt.vatf -q ../data/test_forward_reads.fa -a match --whole-genome
VAT dna -d hp.vatf -q ../data/H_pyloriJ99_Eslice.fasta -a match --whole-genome --match 1 --mismatch -1
VAT view -a match.vatr -o match
vim match
```
## Protein
```console
VAT makevatdb --dbtype prot --in ../data/pfam_ref.fa -d pfam
VAT protein -d pfam.vatf -q ../data/pfam_test.fa -a match -p 4
VAT view -a match.vatr -o match
vim match
```
