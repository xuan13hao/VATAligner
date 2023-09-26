VAT dna -d marine.ref.vatf -q ~/test_data/test_data/marine.sim.singular.fn/marine.sim.singular.fn -a marine.singular.match --for_only
VAT dna -d marine.ref.vatf -q marine.sim.duplex.fn -a marine.duplex.match --for_only --chimera --match 5 --mismatch -4
VAT dna -d Anabas_testudineus.fAnaTes1.2.dna.toplevel.vatf -q Anabas_exons_sim.fa -a Anabas.splice.match --for_only --splice
VAT dna -d marine.ref.vatf -q marine.sim.long.fa -a marine.long.match --for_only --splice --long-read
VAT dna -d human.ref.vatf -q chim.fn -a human.chim.match --whole-genome
VAT protein -d pfam.vatf -q pfam_test.fa -a pfam.match
VAT protein -d pfam.vatf -q pfam_test.fa -a pfam.match --accuracy