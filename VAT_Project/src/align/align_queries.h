

#ifndef ALIGN_QUERIES_H_
#define ALIGN_QUERIES_H_

#include <memory>
#include "../tools/merge_sort.h"
#include "../basic/DiagonalSegment.h"
#include "../search/trace_pt_buffer.h"
#include "../tools/map.h"
#include "align_read.h"
#include "../tools/task_queue.h"
#include "../tools/Queue.h"
#include "../basic/Hits.h"
#include "./Align_fectcher.h"
#include "../search/XdropUngapped.h"

#define MAX_CONTEXT 6
using std::vector;
using std::unique_ptr;

template<typename _val>
struct WorkTarget {
	WorkTarget(size_t block_id, const sequence<const _val>  &seq) :
		block_id(block_id),
		seq(seq),
		filter_score(0),
		outranked(false)
	{}
	bool operator<(const WorkTarget &t) const {
		return filter_score > t.filter_score || (filter_score == t.filter_score && block_id < t.block_id);
	}
	size_t block_id;
	sequence<const _val> seq;
	int filter_score;
	bool outranked;
	// std::array<std::list<Hsp_traits>, MAX_CONTEXT> hsp;
};

struct SeedHit {
	int i, j;
	unsigned frame;
};
struct Output_writer
{
	Output_writer(OutputStreamer* f):
		f_ (f)
	{ }
	void operator()(Text_buffer &buf)
	{
		f_->write(buf.get_begin(), buf.size());
		buf.clear();
	}
private:
	OutputStreamer* const f_;
};
template<typename _val>
void ungapped_stage(const SeedHit *begin, const SeedHit *end, const sequence<const _val> query_seq, size_t block_id) 
{
	array<vector<DiagonalSegment>, MAX_CONTEXT> diagonal_segments;
	WorkTarget<_val> target(block_id, ReferenceSeqs<_val>::get()[block_id]);
	for (const SeedHit *hit = begin; hit < end; ++hit) 
	{
		// const auto d = xdrop_ungapped(query_seq[hit->frame], target.seq, hit->i, hit->j);
		// if(d.score >= 23) diagonal_segments[hit->frame].push_back(d);
	}
	for (unsigned frame = 0; frame < 1; ++frame) 
	{
		// std::stable_sort(diagonal_segments[frame].begin(), diagonal_segments[frame].end(), DiagonalSegment::cmp_diag);
	}
}

template<typename _val, typename _locr, typename _locl, unsigned _d>
void alignQueries(typename Trace_pt_list<_locr,_locl>::iterator begin,
		typename Trace_pt_list<_locr,_locl>::iterator end,
		Output_buffer<_val> &buffer,
		Statistics &st)
{
	array<vector<DiagonalSegment>, MAX_CONTEXT> diagonal_segments;
	typedef Map<typename vector<Hits<_locr,_locl> >::iterator,typename Hits<_locr,_locl>::template Query_id<_d> > Map_t;
	Map_t hits_ (begin, end);
	typename Map_t::Iterator i = hits_.begin();
	const SequenceSet<_val> *ref = ReferenceSeqs<_val>::data_;
	while (i.valid())
	{
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator b = i.begin();
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator e = i.end();
		for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator c = b; c != e; ++c) 
		{
			const _val* subject = ReferenceSeqs<_val>::data_->data(c->subject_);
			const _val* query = QuerySeqs<_val>::data_->data(c->query_);
			const auto d = xdrop_ungapped<_val,_locr,_locl>(query,subject,c->query_,c->subject_);
			if(d.score >= VATParameters::min_ungapped_raw_score) diagonal_segments[c->query_].push_back(d);
			cout<<"subject_ = "<<c->subject_<<", query = "<<c->query_<<", offset = "<<(int)c->seed_offset_<<endl;

		}
		++i;
	}
	while(i.valid() && !exception_state()) 
	{
		align_read<_val,_locr,_locl>(buffer, st, i.begin(), i.end());
		++i;
	}
}

template<typename _val, typename _locr, typename _locl, typename _buffer>
struct Align_context
{
	Align_context(Trace_pt_list<_locr,_locl> &trace_pts, OutputStreamer* output_file):
		trace_pts (trace_pts),
		output_file (output_file),
		writer (output_file),
		queue (VATParameters::threads()*8, writer)
	{ }
	void operator()(unsigned thread_id)
	{
		Statistics st;
		size_t i=0;
		//vector to store
		typename Trace_pt_list<_locr,_locl>::Query_range query_range (trace_pts.get_range());
		_buffer *buffer = 0;
		while(queue.get(i, buffer, query_range) && !exception_state()) 
		{
			try 
			{

				alignQueries<_val,_locr,_locl,1>(query_range.begin, query_range.end, *buffer, st);
				queue.push(i);
			} catch(std::exception &e) {
				exception_state.set(e);
				queue.wake_all();
			}
		}
		statistics += st;
	}
	Trace_pt_list<_locr,_locl> &trace_pts;
	OutputStreamer* output_file;
	Output_writer writer;
	Task_queue<_buffer,Output_writer> queue;
};

template<typename _val, typename _locr, typename _locl>
void alignQueries(const Trace_pt_buffer<_locr,_locl> &trace_pts, OutputStreamer* output_file)
{
	Trace_pt_list<_locr,_locl> v;
	pair<size_t, size_t> query_range;
	for(unsigned bin=0;bin<trace_pts.bins();++bin) 
	{
		// log_stream << "Processing query bin " << bin+1 << '/' << trace_pts.bins() << '\n';
		TimerTools timer ("Loading trace points", false);
		trace_pts.load(v, bin);
		merge_sort(v.begin(), v.end(), VATParameters::threads());
		v.init();

		// Align_fetcher<_locr,_locl>::init(query_range.first, query_range.second, v.begin(), v.end());


		timer.go("Computing alignments");
		if(ref_header.n_blocks > 1) {
			Align_context<_val,_locr,_locl,Temp_output_buffer<_val> > context (v, output_file);
			launch_thread_pool(context, VATParameters::threads());
		} else {
			Align_context<_val,_locr,_locl,Output_buffer<_val> > context (v, output_file);
			launch_thread_pool(context, VATParameters::threads());
		}
	}
	// cout<<"align_queries"<<endl;
}

template<typename _val, typename _locr, typename _locl>
void align(const Trace_pt_buffer<_locr,_locl> &trace_pts, OutputStreamer* output_file)
{
	Consumer* outputfile;
	Trace_pt_list<_locr,_locl> v;
	pair<size_t, size_t> query_range;
	for(unsigned bin=0;bin<trace_pts.bins();++bin) 
	{
		// log_stream << "Processing query bin " << bin+1 << '/' << trace_pts.bins() << '\n';
		TimerTools timer ("Loading trace points", false);
		trace_pts.load(v, bin);
		merge_sort(v.begin(), v.end(), VATParameters::threads());
		v.init();

		// Align_fetcher<_locr,_locl>::init(query_range.first, query_range.second, v.begin(), v.end());

		// OutputSink::instance = unique_ptr<OutputSink>(new OutputSink(query_range.first, outputfile));

		timer.go("Computing alignments");
		if(ref_header.n_blocks > 1) {
			Align_context<_val,_locr,_locl,Temp_output_buffer<_val> > context (v, output_file);
			launch_thread_pool(context, VATParameters::threads());
		} else {
			Align_context<_val,_locr,_locl,Output_buffer<_val> > context (v, output_file);
			launch_thread_pool(context, VATParameters::threads());
		}
	}
	// cout<<"align_queries"<<endl;
}


#endif /* ALIGN_QUERIES_H_ */
