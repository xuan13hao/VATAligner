

#ifndef ALIGN_SEQUENCE_H_
#define ALIGN_SEQUENCE_H_

#include <vector>
#include "../dp/floating_sw.h"

using std::vector;

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
	std::sort(begin, end, hit<_locr,_locl>::cmp_normalized_subject);
	const unsigned q_num (begin->query_);
	const sequence<const _val> query (query_seqs<_val>::get()[q_num]);
	const unsigned frame = q_num % query_contexts();
	const unsigned query_len = query.length();
	padding[frame] = VATParameters::read_padding<_val>(query_len);
	cout<<"align_sequence 1"<<endl;

	const SequenceSet<_val> *ref = ref_seqs<_val>::data_;
	for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator i = begin; i != end; ++i) 
	{
		cout<<"align_sequence 2"<<endl;
		if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) {
			stat.inc(Statistics::DUPLICATES);
			continue;
		}
		local.push_back(local_match<_val> (i->seed_offset_, ref->data(i->subject_)));
		floating_sw(&query[i->seed_offset_],
				local.back(),
				padding[frame],
				score_matrix::get().rawscore(VATParameters::gapped_xdrop),
				VATParameters::gap_open + VATParameters::gap_extend,
				VATParameters::gap_extend,
				transcript_buf,
				Traceback ());
		const int score = local.back().score_;
		std::pair<size_t,size_t> l = ref_seqs<_val>::data_->local_position(i->subject_);
		matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
		anchored_transform(local.back(), l.second, i->seed_offset_);
		stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);

		//local.back().print(query, ref_seqs<_val>::get()[l.first], transcript_buf);

		to_source_space(local.back(), frame, dna_len);
		stat.inc(Statistics::SCORE_TOTAL, local.back().score_);
		stat.inc(Statistics::OUT_HITS);
			

	}

	cout<<"align_sequence 3"<<endl;

}

#endif /* ALIGN_SEQUENCE_H_ */
