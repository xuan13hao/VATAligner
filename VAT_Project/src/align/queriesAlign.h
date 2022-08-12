

#ifndef ALIGN_QUERIES_H_
#define ALIGN_QUERIES_H_
#include <vector>
#include <assert.h>
#include <memory>
#include "../tools/merge_sort.h"
#include "../basic/DiagonalSeeds.h"
#include "../search/trace_pt_buffer.h"
#include "../tools/map.h"
#include "../tools/task_queue.h"
#include "../tools/Queue.h"
#include "../basic/Hits.h"
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
void alignSequence(vector<Segment<_val> > &matches,
		Statistics &stat,
		vector<local_match<_val> > &local,
		unsigned *padding,
		size_t db_letters,
		unsigned dna_len,
		typename Trace_pt_buffer::Vector::iterator &begin,
		typename Trace_pt_buffer::Vector::iterator &end,
		vector<char> &transcript_buf)
{

	std::sort(begin, end, hit::cmp_normalized_subject);
	const unsigned q_num (begin->query_);
	const sequence<const _val> query (QuerySeqs<_val>::get()[q_num]);
	const unsigned frame = q_num % query_contexts();
	const unsigned query_len = query.length();
	padding[frame] = VATParameters::read_padding<_val>(query_len);
	const SequenceSet<_val> *ref = ReferenceSeqs<_val>::data_;

	vector<DiagonalSeeds> diagonalsegment_;
	for(typename Trace_pt_buffer::Vector::iterator i = begin; i != end; ++i) 
	{
		if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) 
		{
			stat.inc(Statistics::DUPLICATES);
			continue;
		}
		std::pair<size_t, size_t> l = ref->local_position(i->subject_);
		const _val* sbj = ref->data(i->subject_);
		const _val* qry = &query[i->seed_offset_];
		// cout<<"i->query_ = "<<(int)i->query_<<", l.second = "<<(int)l.second<<", i->subject = "<<i->subject_<<endl;
		DiagonalSeeds ds = xdrop_ungapped<_val>(qry, sbj,(int)i->seed_offset_,(int)l.second,*i);
		diagonalsegment_.push_back(ds);
	}
	cout<<"diagonalsegment size = "<<diagonalsegment_.size()<<endl;
	for (size_t i = 0; i < diagonalsegment_.size(); i++)
	{
		cout<<"i = "<<diagonalsegment_[i].i<<", j = "<<diagonalsegment_[i].j<<", len = "<<diagonalsegment_[i].len<<", score = "<<diagonalsegment_[i].score<<endl;
	}


	for(typename Trace_pt_buffer::Vector::iterator i = begin; i != end; ++i) 
	{
		if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) 
		{
			stat.inc(Statistics::DUPLICATES);
			continue;
		}

		local.push_back(local_match<_val> (i->seed_offset_, ref->data(i->subject_)));
		floatingSmithWaterman(&query[i->seed_offset_],
				local.back(),
				padding[frame],
				ScoreMatrix::get().rawscore(VATParameters::gapped_xdrop),
				VATParameters::gap_open + VATParameters::gap_extend,
				VATParameters::gap_extend,
				transcript_buf,
				Traceback ());
		const int score = local.back().score_;
		std::pair<size_t,size_t> l = ReferenceSeqs<_val>::data_->local_position(i->subject_);
		matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
		anchored_transform(local.back(), l.second, i->seed_offset_);
		stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);
		// local.back().print(query, ref_seqs<_val>::get()[l.first], transcript_buf);
		to_source_space(local.back(), frame, dna_len);
		stat.inc(Statistics::SCORE_TOTAL, local.back().score_);
		stat.inc(Statistics::OUT_HITS);
	
	}
}

/**
 * @brief Locate read for align sequence 
 * 
 * @tparam _val char DNA Protein
 * @tparam _locr uint64
 * @tparam _locl uint64
 * @param buffer 
 * @param stat statistic 
 * @param begin begin of hit
 * @param end end of hit
 */


template<typename _val>
void alignRead(Output_buffer<_val> &buffer,
		Statistics &stat,
		typename Trace_pt_buffer::Vector::iterator &begin,
		typename Trace_pt_buffer::Vector::iterator &end
		)
{
	static thread_specific_ptr<vector<local_match<_val> > > local_ptr;
	static thread_specific_ptr<vector<Segment<_val> > > matches_ptr;
	static thread_specific_ptr<vector<char> > transcript_ptr;

	Tls<vector<Segment<_val> > > matches (matches_ptr);
	Tls<vector<local_match<_val> > > local (local_ptr);
	Tls<vector<char> > transcript_buf (transcript_ptr);
	local->clear();
	matches->clear();
	transcript_buf->clear();
	assert(end > begin);
	const size_t hit_count = end - begin;
	local->reserve(hit_count);
	const unsigned contexts = query_contexts();
	const unsigned query = begin->query_/contexts;
	const size_t query_len (QuerySeqs<_val>::data_->length(query*contexts));
	const size_t source_query_len = query_len;
	// const size_t source_query_len = query_translated() ? query_seqs<_val>::data_->reverse_translated_len(query*contexts) : query_len;
	const size_t db_letters = ref_header.letters;
	unsigned padding[6];
	typedef Map<typename vector<hit >::iterator,typename hit::template Query_id<1> > Map_t;
	Map_t hits_ (begin, end);
	typename Map_t::Iterator i = hits_.begin();

	while(i.valid()) 
	{

		alignSequence<_val>(*matches, stat, *local, padding, db_letters, source_query_len, i.begin(), i.end(), *transcript_buf);
		++i;
	}


	if(matches->size() == 0)
		return;

	link_segments(*matches);


	std::sort(matches->begin(), matches->end());
	unsigned n_hsp = 0, n_target_seq = 0;
	typename vector<Segment<_val> >::iterator it = matches->begin();

	
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
		// cout<<"frame = "<<it->frame_<<endl;
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

}


template<typename _val,unsigned _d>
void alignQueries(typename Trace_pt_list::iterator begin,
		typename Trace_pt_list::iterator end,
		Output_buffer<_val> &buffer,
		Statistics &st)
{
	// vector<DiagonalSeeds> diagonal_segments;
	typedef Map<typename vector<hit>::iterator,typename hit::template Query_id<_d> > Map_t;
	Map_t hits_ (begin, end);
	typename Map_t::Iterator i = hits_.begin();
	
	const SequenceSet<_val> *ref = ReferenceSeqs<_val>::data_;

	while(i.valid() && !exception_state()) 
	{
		alignRead<_val>(buffer, st, i.begin(), i.end());
		++i;
	}
}

template<typename _val, typename _buffer>
struct Align_context
{
	Align_context(Trace_pt_list &trace_pts, OutputStreamer* output_file):
		trace_pts (trace_pts),
		output_file (output_file),
		writer (output_file),
		queue (VATParameters::threads()*8, writer)
	{ }
	void operator()(unsigned thread_id)
	{
		Statistics st;
		size_t i=0;
		typename Trace_pt_list::Query_range query_range (trace_pts.get_range());
		_buffer *buffer = 0;
		while(queue.get(i, buffer, query_range) && !exception_state()) 
		{
			try 
			{

				alignQueries<_val,1>(query_range.begin, query_range.end, *buffer, st);
				queue.push(i);
			} catch(std::exception &e) {
				exception_state.set(e);
				queue.wake_all();
			}
		}
		statistics += st;
	}
	Trace_pt_list &trace_pts;
	OutputStreamer* output_file;
	Output_writer writer;
	Task_queue<_buffer,Output_writer> queue;
};

template<typename _val, typename _locr, typename _locl>
void alignQueries(Trace_pt_buffer &trace_pts, OutputStreamer* output_file)
{

	// vector<sequence< _val> > subjects;
	// subjects.reserve(ReferenceSeqs<_val>::get().get_length());
	// for (size_t i = 0; i < ReferenceSeqs< _val>::get().get_length(); ++i)
	// 	subjects.push_back(ReferenceSeqs<_val>::get()[i]);
	
	// cout<<"size = "<<subjects.size()<<endl;
	Trace_pt_list v;
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
			Align_context<_val,Temp_output_buffer<_val> > context (v, output_file);
			launch_thread_pool(context, VATParameters::threads());
		} else {
			Align_context<_val,Output_buffer<_val> > context (v, output_file);
			launch_thread_pool(context, VATParameters::threads());
		}
	}
}


#endif /* ALIGN_QUERIES_H_ */
