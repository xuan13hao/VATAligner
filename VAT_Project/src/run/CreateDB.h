
#ifndef MAKE_DB_H_
#define MAKE_DB_H_

#include <iostream>
#include "../basic/VATParameters.h"
#include "../data/Reference.h"
#include "../basic/exceptions.h"
#include "../basic/statistics.h"
#include "../data/LoadSequence.h"
#include "../util/seq_file_format.h"

template<class _val>
void CreateDB(_val)
{
	using std::cout;
	using std::endl;

	cout << "Database File = " << VATParameters::input_ref_file << endl;
	boost::timer::cpu_timer total;
	task_timer timer ("Opening the database file", true);
	Input_stream db_file (VATParameters::input_ref_file, true);
	timer.finish();
	ref_header.block_size = VATParameters::chunk_size;
	size_t chunk = 0;
	Output_stream main(VATParameters::database);
	main.write(&ref_header, 1);

	for(;;++chunk) {
		timer.go("Loading sequences");
		SequenceSet<DNA>* ss;
		size_t n_seq = load_seqs<_val,_val,Single_strand>(db_file, FASTA_format<_val> (), (SequenceSet<_val>**)&ref_seqs<_val>::data_, ref_ids::data_, ss, (size_t)(VATParameters::chunk_size * 1e9));
		if(n_seq == 0)
			break;
		ref_header.letters += ref_seqs<_val>::data_->letters();
		ref_header.sequences += n_seq;
		const bool long_addressing = ref_seqs<_val>::data_->raw_len() > (size_t)std::numeric_limits<uint32_t>::max();
		ref_header.long_addressing = ref_header.long_addressing == true ? true : long_addressing;
		timer.finish();
		ref_seqs<_val>::data_->print_stats();

		timer.go("Building histograms");
		SeedHistogram *hst = new SeedHistogram (*ref_seqs<_val>::data_, _val());

		timer.go("Saving to disk");
		ref_seqs<_val>::data_->save(main);
		ref_ids::get().save(main);
		hst->save(main);

		timer.go("Deallocating sequences");
		delete ref_seqs<_val>::data_;
		delete ref_ids::data_;
		delete ss;
		delete hst;
	}

	timer.finish();
	ref_header.n_blocks = chunk;
	main.seekp(0);
	main.write(&ref_header, 1);
	main.close();

	cout << "Total time = " << boost::timer::format(total.elapsed(), 1, "%ws\n");
}

#endif /* MAKE_DB_H_ */
