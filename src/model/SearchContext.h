#ifndef __SEARCHCONTEXT_H__
#define __SEARCHCONTEXT_H__

#include <iostream>
#include <boost/timer/timer.hpp>
#include "../Database/Reference.h"
#include "../Database/Queries.h"
#include "../Commons/Statistics.h"
#include "../Commons/ShapeParameter.h"
#include "../Output/join_blocks.h"
#include "../Alignment/queriesAlign.h"
#include "../Filter/AlignPartition.h"
#include "../Commons/ContextSet.h"

using std::endl;
using std::cout;
using boost::timer::cpu_timer;
using boost::ptr_vector;


template<typename _val, typename _locr, typename _locq, typename _locl>
class Search_context
{
    public:
	Search_context(unsigned sid, const typename SortedList<_locr>::Type &ref_idx, const typename SortedList<_locq>::Type &query_idx):
		sid (sid),
		ref_idx (ref_idx),
		query_idx (query_idx)
	{ }
	void operator()(unsigned thread_id, unsigned seedp) const
	{
		Statistics stat;
		alignPartition<_val,_locr,_locq,_locl>(seedp,
				stat,
				sid,
				ref_idx.get_partition_cbegin(seedp),
				query_idx.get_partition_cbegin(seedp),
				thread_id);
		statistics += stat;
	}
	const unsigned sid;
	const typename SortedList<_locr>::Type &ref_idx;
	const typename SortedList<_locq>::Type &query_idx;
};

template<typename _val, typename _locr, typename _locq, typename _locl>
void processShapes(unsigned sid,
		cpu_timer &timer_mapping,
		unsigned query_chunk,
		char *query_buffer,
		char *ref_buffer)
{
	using std::vector;

	// partition p (VATConsts::seedp, VATParameters::lowmem);
	::partition p (VATConsts::seedp, VATParameters::lowmem);
	for(unsigned chunk=0;chunk < p.parts; ++chunk) 
	{

		const seedp_range range (p.getMin(chunk), p.getMax(chunk));
		// cout<<"p.getMin(chunk) = "<<p.getMin(chunk)<<", p.getMax(chunk) = "<<p.getMax(chunk)<<endl;
		current_range = range;
		TimerTools timer ("Building reference index", true);
		typename SortedList<_locr>::Type ref_idx (ref_buffer,
				*ReferenceSeqs<_val>::data_,
				ShapeConfigures::instance.get_shape(sid),
				ref_hst.get(VATParameters::index_mode, sid),
				range);
		ReferenceSeqs<_val>::get_nc().template build_masking<_locr>(sid, range, ref_idx);

		timer.go("Building query index");
		timer_mapping.resume();
		typename SortedList<_locq>::Type query_idx (query_buffer,
				*QuerySeqs<_val>::data_,
				ShapeConfigures::instance.get_shape(sid),
				query_hst->get(VATParameters::index_mode, sid),
				range);
		timer.finish();

		timer.go("Searching seeds");
		Search_context<_val,_locr,_locq,_locl> context (sid, ref_idx, query_idx);

		launch_scheduled_thread_pool(context, VATConsts::seedp, VATParameters::thread());

	}
	timer_mapping.stop();
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void ProcessRefsChunks(Database_file<_val> &db_file,
		cpu_timer &timer_mapping,
		cpu_timer &total_timer,
		unsigned query_chunk,
		pair<size_t,size_t> query_len_bounds,
		char *query_buffer,
		VATOutput &master_out,
		vector<TempFile> &tmp_file)
{
	TimerTools timer ("Loading reference sequences", true);
	ReferenceSeqs<_val>::data_ = new Masked_sequence_set<_val> (db_file);
	ReferenceIds::data_ = new AlphabetSet<char,0> (db_file);
	ref_hst.load(db_file);

	setup_search_params<_val>(query_len_bounds, ReferenceSeqs<_val>::data_->letters());

	ref_map.init(ReferenceSeqs<_val>::get().get_length());
	timer.go("Allocating buffers");
	char *ref_buffer = SortedList<_locr>::Type::alloc_buffer(ref_hst);

	timer.go("Initializing temporary storage");
	timer_mapping.resume();
	Trace_pt_buffer<_locr,_locl>::instance = new Trace_pt_buffer<_locr,_locl> (QuerySeqs<_val>::data_->get_length()/query_contexts(),
			VATParameters::tmpdir,
			VATParameters::mem_buffered());
	timer.finish();
	timer_mapping.stop();

	for(unsigned i=0;i<ShapeConfigures::instance.count();++i)
		processShapes<_val,_locr,_locq,_locl>(i, timer_mapping, query_chunk, query_buffer, ref_buffer);

	timer.go("Closing temporary storage");
	Trace_pt_buffer<_locr,_locl>::instance->close();
	exception_state.sync();
	timer.go("Deallocating buffers");
	delete[] ref_buffer;

	timer_mapping.resume();
	OutputStreamer* out;
	if(ref_header.n_blocks > 1) {
		timer.go ("Opening temporary output file");
		tmp_file.push_back(TempFile ());
		out = new OutputStreamer (tmp_file.back());
	} else
		out = &master_out.stream();
	timer.go("Generating seeds...");
	alignQueries<_val,_locr,_locl>(*Trace_pt_buffer<_locr,_locl>::instance, out);
	delete Trace_pt_buffer<_locr,_locl>::instance;

	if(ref_header.n_blocks > 1) {
		timer.go("Closing temporary output file");
		out->close();
		delete out;
	}
	timer_mapping.stop();

	timer.go("Deallocating reference");
	delete ReferenceSeqs<_val>::data_;
	delete ReferenceIds::data_;
	timer.finish();
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void ProcessQueryChunks(Database_file<_val> &db_file,
		cpu_timer &timer_mapping,
		cpu_timer &total_timer,
		unsigned query_chunk,
		pair<size_t,size_t> query_len_bounds,
		VATOutput &master_out)
{
	TimerTools timer ("Allocating buffers", true);
	char *query_buffer = SortedList<_locq>::Type::alloc_buffer(*query_hst);
	vector<TempFile> tmp_file;
	timer.finish();
	db_file.rewind();
	for(current_ref_block=0;current_ref_block<ref_header.n_blocks;++current_ref_block)
		ProcessRefsChunks<_val,_locr,_locq,_locl>(db_file, timer_mapping, total_timer, query_chunk, query_len_bounds, query_buffer, master_out, tmp_file);

	timer.go("Deallocating buffers");
	timer_mapping.resume();
	delete[] query_buffer;

	if(ref_header.n_blocks > 1) 
	{
		timer.go("Joining output blocks");
		join_blocks<_val>(ref_header.n_blocks, master_out, tmp_file);
	}

	timer.go("Deallocating queries");
	delete QuerySeqs<_val>::data_;
	delete query_ids::data_;
	delete query_source_seqs::data_;
	timer_mapping.stop();
}




#endif // __SEARCHCONTEXT_H__