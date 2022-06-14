

#ifndef ALIGN_READ_H_
#define ALIGN_READ_H_

#include <vector>
#include <assert.h>
#include "../tools/async_buffer.h"
#include "../basic/Hits.h"
#include "../basic/statistics.h"
#include "../search/align_ungapped.h"
#include "align_sequence.h"
#include "../tools/text_buffer.h"
#include "../output/output_buffer.h"
#include "link_segments.h"

using std::vector;

template<typename _val, typename _locr, typename _locl>
void align_read(Output_buffer<_val> &buffer,
		Statistics &stat,
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &begin,
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &end)
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


	typedef Map<typename vector<Hits<_locr,_locl> >::iterator,typename Hits<_locr,_locl>::template Query_id<1> > Map_t;
	Map_t hits (begin, end);

	typename Map_t::Iterator i = hits.begin();
	while(i.valid()) {

		align_sequence<_val,_locr,_locl>(*matches, stat, *local, padding, db_letters, source_query_len, i.begin(), i.end(), *transcript_buf);
		++i;
	}


	if(matches->size() == 0)
		return;

	link_segments(*matches);


	std::sort(matches->begin(), matches->end());
	unsigned n_hsp = 0, n_target_seq = 0;
	typename vector<Segment<_val> >::iterator it = matches->begin();
	// const int min_raw_score = 10;
	const int min_raw_score = ScoreMatrix::get().rawscore(VATParameters::min_bit_score == 0
			? ScoreMatrix::get().bitscore(VATParameters::max_evalue, ref_header.letters, query_len) : VATParameters::min_bit_score);
	cout<<"min_raw_score = "<<min_raw_score<<endl;
	const int top_score = matches->operator[](0).score_;


	while(it < matches->end()) {
		const bool same_subject = it != matches->begin() && (it-1)->subject_id_ == it->subject_id_;
		if(!same_subject && it->score_ < min_raw_score)
			break;
		if(!same_subject && !VATParameters::output_range(n_target_seq, it->score_, top_score))
			break;
		if(same_subject && (it-1)->score_ == it->score_) {
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

}

#endif /* ALIGN_READ_H_ */
