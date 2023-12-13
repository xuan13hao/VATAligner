make
./VAT makevatdb --dbtype nucl --in /home/xuan/test_data/test_data/marine.ref.fn -d marine.ref -b 1
/usr/bin/time -v ./VAT dna -d marine.ref -q /home/xuan/test_data/test_data/marine.sim.singular.fn/marine.sim.singular.fn -a marine.ref -p 4 -t /dev/shm --for_only -o match --id 14
perl ../scripts/Benchmark_DIAMOND_Singular.pl /home/xuan/test_data/test_data/marine.sim.singular.fn/marine.sim.singular.fn  match
