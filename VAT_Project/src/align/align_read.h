

#ifndef ALIGN_READ_H_
#define ALIGN_READ_H_

#include <vector>
#include <assert.h>
#include "../tools/async_buffer.h"
#include "../basic/Hits.h"
#include "../basic/Statistics.h"
#include "../search/XdropUngapped.h"
#include "../tools/text_buffer.h"
#include "../output/output_buffer.h"
#include "link_segments.h"
#include "../dp/floating_sw.h"

using std::vector;
struct SeedHit {
	int i, j;
	unsigned frame;
};
/**
 * @brief Align Sequence
 * 
 * @tparam _val char DNA or Protein
 * @tparam _locr uint_64
 * @tparam _locl uint_64
 * @param matches vector for segement
 * @param stat statistic 
 * @param local local match
 * @param padding 
 * @param db_letters 
 * @param dna_len 
 * @param begin begin of Hits
 * @param end end of Hits
 * @param transcript_buf Edit_transcript
 */


template<typename _val, typename _locr, typename _locl>
void align_sequence(vector<Segment<_val> > &matches,
		Statistics &stat,
		vector<local_match<_val> > &local,
		unsigned *padding,
		size_t db_letters,
		unsigned dna_len,
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &begin,
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &end,
		vector<char> &transcript_buf)
{

	std::sort(begin, end, Hits<_locr,_locl>::cmp_normalized_subject);
	const unsigned q_num (begin->query_);
	const sequence<const _val> query (QuerySeqs<_val>::get()[q_num]);
	const unsigned frame = q_num % query_contexts();
	const unsigned query_len = query.length();
	padding[frame] = VATParameters::read_padding<_val>(query_len);
	const SequenceSet<_val> *ref = ReferenceSeqs<_val>::data_;
	/*
	vector<SeedHit>* seed_hits;
	for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator c = begin; c != end; ++c) 
	{
		std::pair<size_t, size_t> l = ref->local_position(c->subject_);
		// if (l.first != target) 
		// {
		// 	hits.next();
		// 	target = l.first;
		// 	target_block_ids.push_back(target);
		// }
		cout<<" seed_offset_= "<<(int)c->seed_offset_<<", s = "<<(int)l.second<<", q = "<<(int)c->query_<<endl;
		seed_hits.push_back({ (int)c->seed_offset_, (int)l.second, (int)c->query_ });
		// seed_hits.push_back({ (int64_t)c->query_, (int64_t)c->subject_, (int64_t)c->seed_offset_ });

	}
	cout<<"hit size = "<<seed_hits.size()<<endl;
	for (size_t i = 0; i < seed_hits.size(); i++)
	{
		xdrop_ungapped();
	}
	*/
	vector<DiagonalSegment> diagonalsegment_;

	for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator i = begin; i != end; ++i) 
	{
		// if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) 
		// {
		// 	stat.inc(Statistics::DUPLICATES);
		// 	continue;
		// }
		std::pair<size_t, size_t> l = ref->local_position(i->subject_);
		const _val* sbj = ref->data(i->subject_);
		const _val* qry = &query[i->seed_offset_];
		// cout<<"i->query_ = "<<(int)i->query_<<", l.second = "<<(int)l.second<<", i->subject = "<<i->subject_<<endl;
		DiagonalSegment ds = xdrop_ungapped<_val, _locr,_locl>(qry, sbj,(int)i->seed_offset_,(int)l.second);
		diagonalsegment_.push_back(ds);
	}
	cout<<"diagonalsegment size = "<<diagonalsegment_.size()<<endl;

	for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator i = begin; i != end; ++i) 
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


template<typename _val, typename _locr, typename _locl>
void align_read(Output_buffer<_val> &buffer,
		Statistics &stat,
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &begin,
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &end
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
	typedef Map<typename vector<Hits<_locr,_locl> >::iterator,typename Hits<_locr,_locl>::template Query_id<1> > Map_t;
	Map_t hits_ (begin, end);
	typename Map_t::Iterator i = hits_.begin();

	while(i.valid()) 
	{

		align_sequence<_val,_locr,_locl>(*matches, stat, *local, padding, db_letters, source_query_len, i.begin(), i.end(), *transcript_buf);
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

#endif /* ALIGN_READ_H_ */
