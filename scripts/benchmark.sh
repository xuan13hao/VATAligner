/usr/bin/time -v VAT dna -d marine.ref.vatf -q ~/test_data/test_data/marine.sim.singular.fn/marine.sim.singular.fn -a marine.singular.match --for_only
/usr/bin/time -v VAT dna -d marine.ref.vatf -q marine.sim.duplex.fn -a marine.duplex.match --for_only --chimera --match 5 --mismatch -4
/usr/bin/time -v VAT dna -d Anabas_testudineus.fAnaTes1.2.dna.toplevel.vatf -q Anabas_exons_sim.fa -a Anabas.splice.match --for_only --splice
/usr/bin/time -v VAT dna -d marine.ref.vatf -q marine.sim.long.fa -a marine.long.match --for_only --long-read
/usr/bin/time -v VAT dna -d human.ref.vatf -q chim.fn -a human.chim.match --whole-genome
/usr/bin/time -v VAT protein -d pfam.vatf -q pfam_test.fa -a pfam.match
/usr/bin/time -v VAT protein -d pfam.vatf -q pfam_test.fa -a pfam.match --accuracy