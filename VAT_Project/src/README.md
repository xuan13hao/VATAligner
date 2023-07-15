# Introduction

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

'''R
perfect 8
partial_multi   0
both_multi      0
partial_wrong   66
both_wrong      72844
partial_miss    926930
both_miss       152
>>>Strand-level statistics:
correct 927012
unique  997898

>>>Read-level statistics:
perfect 795095
partial_multi   46942
both_multi      295
partial_wrong   31895
both_wrong      111865
partial_miss    13844
both_miss       64
>>>Strand-level statistics:
correct 1682871
unique  1815362

'''

