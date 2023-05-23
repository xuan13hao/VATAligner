
#ifndef MAKE_DB_H_
#define MAKE_DB_H_

#include <iostream>
#include "../basic/VATParameters.h"
#include "../data/Reference.h"
#include "../basic/exceptions.h"
#include "../basic/Statistics.h"
#include "../data/LoadSequence.h"
#include "../tools/seq_file_format.h"

template<class _val>
void CreateDB(_val)
{
	using std::cout;
	using std::endl;

	cout << "Database File = " << VATParameters::input_ref_file << endl;
	boost::timer::cpu_timer total;
	TimerTools timer ("Opening the database file", true);
	Input_stream db_file (VATParameters::input_ref_file, true);
	timer.finish();
	ref_header.block_size = VATParameters::chunk_size;
	int chunk = 0;
	OutputStreamer main(VATParameters::database);
	main.write(&ref_header, 1);

	for(;;++chunk) 
	{
		timer.go("Loading sequences");
		SequenceSet<DNA>* ss;
		int n_seq = 0;
		if(VATParameters::db_type == "nucl")
		{
			n_seq = ReadingSeqs<_val,_val,Double_strand>(db_file, FASTA_format<_val> (), (SequenceSet<_val>**)&ReferenceSeqs<_val>::data_, ReferenceIds::data_, ss, (size_t)(VATParameters::chunk_size * 1e9));
		}else
		{
			n_seq = ReadingSeqs<_val,_val,Single_strand>(db_file, FASTA_format<_val> (), (SequenceSet<_val>**)&ReferenceSeqs<_val>::data_, ReferenceIds::data_, ss, (size_t)(VATParameters::chunk_size * 1e9));
		}
		if(n_seq == 0)
			break;
		ref_header.letters += ReferenceSeqs<_val>::data_->letters();
		ref_header.sequences += n_seq;
		const bool long_addressing = ReferenceSeqs<_val>::data_->raw_len() > (size_t)std::numeric_limits<uint32_t>::max();
		ref_header.long_addressing = ref_header.long_addressing == true ? true : long_addressing;
		timer.finish();
		ReferenceSeqs<_val>::data_->print_stats();

		timer.go("Building SeedHistograms");
		SeedHistogram *hst = new SeedHistogram (*ReferenceSeqs<_val>::data_, _val());

		timer.go("Saving to disk");
		ReferenceSeqs<_val>::data_->save(main);
		ReferenceIds::get().save(main);
		hst->save(main);

		timer.go("Deallocating sequences");
		delete ReferenceSeqs<_val>::data_;
		delete ReferenceIds::data_;
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

/*
void CreateDNA_DB()
{
	using std::cout;
	using std::endl;

	cout << "Database File = " << VATParameters::input_ref_file << endl;
	boost::timer::cpu_timer total;
	TimerTools timer ("Opening the database file", true);
	Input_stream db_file (VATParameters::input_ref_file, true);
	timer.finish();
	ref_header.block_size = VATParameters::chunk_size;
	int chunk = 0;
	OutputStreamer main(VATParameters::database);
	main.write(&ref_header, 1);
	cout<<"DNA DB"<<endl;
	for(;;++chunk) 
	{
		timer.go("Loading sequences");
		SequenceSet<DNA>* ss;
		int n_seq = 0;
		n_seq = ReadingDNASeqs(db_file, FASTA_format<DNA> (), (SequenceSet<DNA>**)&ReferenceSeqs<DNA>::data_, ReferenceIds::data_, ss, (size_t)(VATParameters::chunk_size * 1e9));
		if(n_seq == 0)
			break;
		ref_header.letters += ReferenceSeqs<DNA>::data_->letters();
		ref_header.sequences += n_seq;
		const bool long_addressing = ReferenceSeqs<DNA>::data_->raw_len() > (size_t)std::numeric_limits<uint32_t>::max();
		ref_header.long_addressing = ref_header.long_addressing == true ? true : long_addressing;
		timer.finish();
		ReferenceSeqs<DNA>::data_->print_stats();

		timer.go("Building SeedHistograms");
		SeedHistogram *hst = new SeedHistogram (*ReferenceSeqs<DNA>::data_, DNA());

		timer.go("Saving to disk");
		ReferenceSeqs<DNA>::data_->save(main);
		ReferenceIds::get().save(main);
		hst->save(main);

		timer.go("Deallocating sequences");
		delete ReferenceSeqs<DNA>::data_;
		delete ReferenceIds::data_;
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
*/

#endif /* MAKE_DB_H_ */
