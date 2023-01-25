# VAT_versions

./VAT makevatdb --dbtype nucl --in ../data/MT.fa -d mt
./VAT makevatdb --dbtype nucl -i ../data/MT.fa -d mt
./VAT dna -d marine.vatf -q ../data/marine.sim.duplex.fn -a match
./VAT view -a match.vatr -o match
./VAT view -a match.vatr -o match -f sam
vim match


./VAT protein -d Pfam.vatf -q /media/xuan/VAT/Pfam.sim.singular.fa -a match

q_num = 0, subject = 24,seed offset = 23
q_num = 1, subject = 85,seed offset = 32
q_num = 0, subject = 9,seed offset = 8
q_num = 0, subject = 15,seed offset = 14
q_num = 1, subject = 81,seed offset = 28
q_num = 0, subject = 28,seed offset = 27
q_num = 0, subject = 27,seed offset = 26
q_num = 0, subject = 8,seed offset = 7
q_num = 0, subject = 14,seed offset = 13
q_num = 1, subject = 78,seed offset = 25
q_num = 1, subject = 79,seed offset = 26
q_num = 1, subject = 86,seed offset = 33
q_num = 0, subject = 31,seed offset = 30
q_num = 0, subject = 29,seed offset = 28
q_num = 0, subject = 34,seed offset = 33
q_num = 0, subject = 36,seed offset = 35
q_num = 1, subject = 87,seed offset = 34
q_num = 1, subject = 88,seed offset = 35

seq1    t2      100.0   68      0       0       1       68      1       68      1.6e-29 108.6
seq2    t2      100.0   68      0       0       1       68      53      120     1.6e-29 108.6
