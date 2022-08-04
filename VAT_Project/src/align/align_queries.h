

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



#include <vector>
#include <assert.h>
#include "../tools/async_buffer.h"

#include "../basic/Statistics.h"
#include "../search/XdropUngapped.h"
#include "../tools/text_buffer.h"
#include "../output/output_buffer.h"
#include "link_segments.h"
#include "../dp/floating_sw.h"


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

struct SHit 
{
	int	query_;
	int	subject_;//reference 
	int	seed_offset_;//seed offset

	int64_t global_diagonal() const
	{ 
		return (int64_t)subject_ - (int64_t)seed_offset_; 
	}
	static bool cmp_normalized_subject(const SHit &lhs, const SHit &rhs)
	{
		const uint64_t x = (uint64_t)lhs.subject_ + (uint64_t)rhs.seed_offset_, y = (uint64_t)rhs.subject_ + (uint64_t)lhs.seed_offset_;
		return x < y || (x == y && lhs.seed_offset_ < rhs.seed_offset_);
	}
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
	array<vector<DiagonalSeeds>, MAX_CONTEXT> diagonal_segments;
	WorkTarget<_val> target(block_id, ReferenceSeqs<_val>::get()[block_id]);
	for (const SeedHit *hit = begin; hit < end; ++hit) 
	{
		cout<<"s = "<<hit->i<<", q = "<<hit->j<<", frame = "<<hit->frame<<endl;
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
	vector<DiagonalSeeds> diagonal_segments;
	typedef Map<typename vector<Hits<_locr,_locl> >::iterator,typename Hits<_locr,_locl>::template Query_id<_d> > Map_t;
	Map_t hits_ (begin, end);
	typename Map_t::Iterator i = hits_.begin();
	
	const SequenceSet<_val> *ref = ReferenceSeqs<_val>::data_;

	while(i.valid() && !exception_state()) 
	{
		align_read<_val,_locr,_locl>(buffer, st, i.begin(), i.end());
		++i;
	}

		/*
	vector<SHit> seed_hits;
	while (i.valid())
	{
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator b = i.begin();
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator e = i.end();
		for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator c = b; c != e; ++c) 
		{
			// std::pair<size_t, size_t> l = ReferenceSeqs<_val>::data_->local_position(c->subject_);
			// if (l.first != target) 
			// {
			// 	hits.next();
			// 	target = l.first;
			// 	target_block_ids.push_back(target);
			// }
			seed_hits.push_back({ (int64_t)c->query_, (int64_t)c->subject_, (int64_t)c->seed_offset_ });
				// const _val* subject = ReferenceSeqs<_val>::data_->data(c->subject_);
			// const _val* query = QuerySeqs<_val>::data_->data(c->query_);
			// const _val q = query[c->seed_offset_];
			// const auto d = xdrop_ungapped<_val,_locr,_locl>(&query[c->seed_offset_],ref->data(c->subject_),c->query_,c->subject_);
			// cout<<"score = "<<d.score<<endl;
			// if(d.score >= VATParameters::min_ungapped_raw_score) diagonal_segments.push_back(d);
			// cout<<"subject_ = "<<c->subject_<<", query = "<<c->query_<<", offset = "<<(int)c->seed_offset_<<endl;
		}
		++i;
	}
	static thread_specific_ptr<vector<local_match<_val> > > local_pt;
	static thread_specific_ptr<vector<Segment<_val> > > matches_pt;
	static thread_specific_ptr<vector<char> > transcript_pt;

	Tls<vector<Segment<_val> > > matche (matches_pt);
	Tls<vector<local_match<_val> > > locals (local_pt);
	Tls<vector<char> > transcript_bufs (transcript_pt);
	locals->clear();
	matche->clear();
	transcript_bufs->clear();

	locals->reserve(seed_hits.size());
	


	const unsigned contexts = query_contexts();
	vector<Segment<_val> >  matches;
	vector<local_match<_val> > local;
	vector<char>  transcript_bufs ;
	std::sort(seed_hits.begin(), seed_hits.end(), SHit::cmp_normalized_subject);
	for (vector<SHit>::iterator i = seed_hits.begin(); i != seed_hits.end(); ++i)
	{
		// cout<<"subject_ = "<<i->subject_<<", query = "<<i->query_<<", offset = "<<(int)i->seed_offset_<<endl;
		const unsigned q = i->query_/contexts;
		const size_t query_len (QuerySeqs<_val>::data_->length(q));
		const size_t source_query_len = query_len;
		const size_t db_letters = ref_header.letters;

		unsigned padding[6];
		const unsigned q_num (i->query_);
		const sequence<const _val> query (QuerySeqs<_val>::get()[q_num]);

		const unsigned frame = q_num % query_contexts();
		const unsigned q_len = query.length();
		padding[frame] = VATParameters::read_padding<_val>(q_len);

		if(i != seed_hits.begin()&& (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) 
		{
			// stat.inc(Statistics::DUPLICATES);
			continue;
		}
	
		local.push_back(local_match<_val> (i->seed_offset_, ref->data(i->subject_)));
		floatingSmithWaterman(&query[i->seed_offset_],
				local.back(),
				padding[frame],
				ScoreMatrix::get().rawscore(VATParameters::gapped_xdrop),
				VATParameters::gap_open + VATParameters::gap_extend,
				VATParameters::gap_extend,
				transcript_bufs,
				Traceback ());
		const int score = local.back().score_;
		std::pair<size_t,size_t> l = ReferenceSeqs<_val>::data_->local_position(i->subject_);
		cout<<"frame = "<<frame<<endl;
		matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
		anchored_transform(local.back(), l.second, i->seed_offset_);
		// stat.inc(Statistics::ALIGNED_QLEN, locals.back().query_len_);
		// local.back().print(query, ref_seqs<_val>::get()[l.first], transcript_buf);
		to_source_space(local.back(), frame, source_query_len);
		// stat.inc(Statistics::SCORE_TOTAL, locals.back().score_);
		// stat.inc(Statistics::OUT_HITS);

	}
	cout<<"match size = "<<matches.size()<<endl;

	if(matches.size() == 0)
		return;

	link_segments(matches);


	std::sort(matches.begin(), matches.end());
	unsigned n_hsp = 0, n_target_seq = 0;
	typename vector<Segment<_val> >::iterator it = matches.begin();
	while (it < matches.end())
	{
		cout<<it->traceback_->len_<<"\t"<<it->traceback_->query_begin_<<"\t"<<it->traceback_->subject_begin_<<"\t"<<it->frame_<<endl;;
		++it;
	}
	
	


	const int min_raw_score = ScoreMatrix::get().rawscore(VATParameters::min_bit_score == 0
			? ScoreMatrix::get().bitscore(VATParameters::max_evalue, ref_header.letters, query_len) : VATParameters::min_bit_score);
		
	const int top_score = matches->operator[](0).score_;


	while(it < matches->end()) 
	{
		const bool same_subject = it != matches->begin() && (it-1)->subject_id_ == it->subject_id_;
		if(!same_subject && it->score_ < min_raw_score)
			break;
		if(!same_subject && !VATParameters::output_range(n_target_seq, it->score_, top_score))
			break;
		if(same_subject && (it-1)->score_ == it->score_) 
		{
			++it;
			continue;
		}
		if(static_cast<double>(it->traceback_->identities_)*100/it->traceback_->len_ < VATParameters::min_id) {
			++it;
			continue;
		}
		if(same_subject && VATParameters::single_domain) {
			++it;
			continue;
		}

		if(n_hsp == 0)
			buffer.write_query_record(query);

		buffer.print_match(*it, source_query_len, QuerySeqs<_val>::get()[query*contexts + it->frame_], query, *transcript_buf);

		++n_hsp;
		if(!same_subject)
			++n_target_seq;
		if(VATParameters::alignment_traceback && it->traceback_->gap_openings_ > 0)
			stat.inc(Statistics::GAPPED);
		++it;
	}

	if(n_hsp > 0)
		buffer.finish_query_record();

	stat.inc(Statistics::OUT_MATCHES, matches->size());
	if(ref_header.n_blocks == 1) {
		stat.inc(Statistics::MATCHES, n_hsp);
		if(n_hsp > 0)
			stat.inc(Statistics::ALIGNED);
	}
*/	
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
		TimerTools timer ("Loading trace points", false);
		trace_pts.load(v, bin);
		merge_sort(v.begin(), v.end(), VATParameters::threads());
		v.init();
		timer.go("Computing alignments");
		if(ref_header.n_blocks > 1) 
		{
			Align_context<_val,_locr,_locl,Temp_output_buffer<_val> > context (v, output_file);
			launch_thread_pool(context, VATParameters::threads());
		} else {
			Align_context<_val,_locr,_locl,Output_buffer<_val> > context (v, output_file);
			launch_thread_pool(context, VATParameters::threads());
		}
	}
}
/*
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
*/

#endif /* ALIGN_QUERIES_H_ */
