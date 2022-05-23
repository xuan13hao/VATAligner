#ifndef __DNAALIGNFLOW_H__
#define __DNAALIGNFLOW_H__



#include <iostream>
#include <boost/timer/timer.hpp>
#include "../data/Reference.h"
#include "../data/Queries.h"
#include "../basic/statistics.h"
#include "../basic/ShapeParameter.h"
#include "../output/join_blocks.h"
#include "../align/align_queries.h"
#include "../search/align_range.h"
#include "../basic/ContextSet.h"
#include "SearchContext.h"

using std::endl;
using std::cout;
using boost::timer::cpu_timer;
using boost::ptr_vector;


template<typename _val, typename _locr>
void DNAMasterThread(Database_file<_val> &db_file, cpu_timer &timer_mapping, cpu_timer &total_timer)
{
	ShapeConfigures::instance = ShapeConfigures (VATParameters::index_mode, _val ());
	task_timer timer ("Opening the input file", true);
	timer_mapping.resume();
	const Sequence_file_format<DNA> *format_n (guess_format<DNA>(VATParameters::query_file));
	Input_stream query_file (VATParameters::query_file, true);
	current_query_chunk=0;
	timer.go("Opening the output file");
	VATOutput master_out;
	timer_mapping.stop();
	timer.finish();
	
	for(;;++current_query_chunk) {
		task_timer timer ("Loading query sequences", true);
		timer_mapping.resume();
		size_t n_query_seqs;

		n_query_seqs = load_seqs<_val,_val,Single_strand>(query_file, *format_n, &query_seqs<_val>::data_, query_ids::data_, query_source_seqs::data_, (size_t)(VATParameters::chunk_size * 1e9));

		if(n_query_seqs == 0)
			break;
		timer.finish();
		query_seqs<_val>::data_->print_stats();

		timer.go("Building query histograms");
		query_hst = auto_ptr<SeedHistogram> (new SeedHistogram (*query_seqs<_val>::data_, _val()));
		const pair<size_t,size_t> query_len_bounds = query_seqs<_val>::data_->len_bounds(ShapeConfigures::get().get_shape(0).length_);
		timer_mapping.stop();
		timer.finish();
		const bool long_addressing_query = query_seqs<_val>::data_->raw_len() > (size_t)std::numeric_limits<uint32_t>::max();

		if(query_len_bounds.second <= (size_t)std::numeric_limits<uint8_t>::max()) {
			if(long_addressing_query)
				run_query_chunk<_val,_locr,uint64_t,uint8_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds, master_out);
			else
				run_query_chunk<_val,_locr,uint32_t,uint8_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds, master_out);
		} else if(query_len_bounds.second <= (size_t)std::numeric_limits<uint16_t>::max()) {
			if(long_addressing_query)
				run_query_chunk<_val,_locr,uint64_t,uint16_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds, master_out);
			else
				run_query_chunk<_val,_locr,uint32_t,uint16_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds, master_out);
		} else {
			if(long_addressing_query)
				run_query_chunk<_val,_locr,uint64_t,uint32_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds, master_out);
			else
				run_query_chunk<_val,_locr,uint32_t,uint32_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds, master_out);
		}
	}

	timer.go("Closing the output file");
	timer_mapping.resume();
	master_out.finish();
	timer_mapping.stop();

	timer.go("Closing the database file");
	db_file.close();

	timer.finish();
	cout << "Total time = " << boost::timer::format(total_timer.elapsed(), 1, "%ws\n");
	cout << "Mapping time = " << boost::timer::format(timer_mapping.elapsed(), 1, "%ws\n");
	statistics.print();
}

template<typename _val>
void DNAMasterThread()
{
	cpu_timer timer2, timer_mapping;
	timer_mapping.stop();

	if(!check_dir(VATParameters::tmpdir))
		throw std::runtime_error("Temporary directory " + VATParameters::tmpdir + " does not exist or is not a directory. Please use option -t to specify a different directory.");

	task_timer timer ("Opening the database", 1);
	Database_file<_val> db_file;
	timer.finish();
	VATParameters::set_options<_val>(ref_header.block_size);
	verbose_stream << "Reference = " << VATParameters::database << endl;
	verbose_stream << "Sequences = " << ref_header.sequences << endl;
	verbose_stream << "Letters = " << ref_header.letters << endl;
	verbose_stream << "Block size = " << (size_t)(ref_header.block_size * 1e9) << endl;

	if(ref_header.long_addressing)
		DNAMasterThread<_val,uint64_t>(db_file, timer_mapping, timer2);
	else
		DNAMasterThread<_val,uint32_t>(db_file, timer_mapping, timer2);
}


#endif // __DNAALIGNFLOW_H__