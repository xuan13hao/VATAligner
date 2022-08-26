#ifndef __ALIGNREADS_H__
#define __ALIGNREADS_H__
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

template<typename _val>
vector<DiagonalSeeds> getSegments(
		typename Trace_pt_buffer::Vector::iterator &begin,
		typename Trace_pt_buffer::Vector::iterator &end,
)
{
	string qry,sbj;
	std::sort(begin, end, hit::cmp_normalized_subject);
	const unsigned q_num (begin->query_);
	const sequence<const _val> query (QuerySeqs<_val>::get()[q_num]);
	const unsigned frame = q_num % query_contexts();
	const unsigned query_len = query.length();
	const SequenceSet<_val> *ref = ReferenceSeqs<_val>::data_;
	vector<DiagonalSeeds> diagonalsegment_;
	for(typename Trace_pt_buffer::Vector::iterator i = begin; i != end; ++i) 
	{
		if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) 
		{
			// stat.inc(Statistics::DUPLICATES);
			continue;
		}
		std::pair<size_t, size_t> l = ref->local_position(i->subject_);
		const _val* sbj = ref->data(i->subject_);
		const _val* qry = &query[i->seed_offset_];
		
		// cout<<"i->query_ = "<<(int)i->query_<<", l.second = "<<(int)l.second<<", i->subject = "<<i->subject_<<endl;
		DiagonalSeeds ds = xdrop_ungapped<_val>(qry, sbj,(int)i->seed_offset_,(int)l.second,*i);
		// const _val* sbj1 = ref->data(ds.i+29);
		// const _val* qry1 = &query[ds.j+28];
		// cout<<"subject = "<<AlphabetAttributes<_val>::ALPHABET[*sbj]<<", query = "<<AlphabetAttributes<_val>::ALPHABET[*qry]<<endl;
		cout<<"i = "<<ds.i<<", j = "<<ds.j<<",len= "<<ds.len<<endl;
		diagonalsegment_.push_back(ds);
	}
    return diagonalsegment_;
}

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

		getSegments<_val>(stat, padding, db_letters, source_query_len, i.begin(), i.end());
		++i;
	}

}
/*
void printMatches(vector<Segment<_val> >& matches,vector<char>&  transcript_buf,int query_len,int query)
{    
    const size_t source_query_len = query_len;
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
		buffer.print_match(*it, source_query_len, QuerySeqs<_val>::get()[query + it->frame_], query, *transcript_buf);

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
*/
#endif // __ALIGNREADS_H__