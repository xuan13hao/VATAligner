

#ifndef ALIGN_SEQUENCE_H_
#define ALIGN_SEQUENCE_H_

#include <vector>
#include "../dp/floating_sw.h"
#include "./pairChimera.h"
#include "./ungappedSeed.h"
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
	std::sort(begin, end, Hits<_locr,_locl>::cmp_normalized_subject);
	const unsigned q_num (begin->query_);
	const sequence<const _val> query (QuerySeqs<_val>::get()[q_num]);
	const unsigned frame = q_num % query_contexts();
	const unsigned query_len = query.length();
	padding[frame] = VATParameters::read_padding<_val>(query_len);

	const SequenceSet<_val> *ref = ReferenceSeqs<_val>::data_;
	
	//get Seeds by ungapped
	vector<DiagonalSeeds<_locr,_locl> > diagonalsegment_;
	for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator i = begin; i != end; ++i)
	{
		if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) 
		{
			stat.inc(Statistics::DUPLICATES);
			continue;
		}
		std::pair<size_t, size_t> l = ref->local_position(i->subject_);
		const _val* sbj = ref->data(i->subject_);
		const _val* qry = &query[i->seed_offset_];
		DiagonalSeeds<_locr,_locl> ds = ungappedSeeds<_val, _locr,_locl> (qry, sbj,(int)i->seed_offset_,(int)l.second,*i);
		if (ds.len > 15)
		{
			diagonalsegment_.push_back(ds);
		}
	}
	
	if(VATParameters::chimera)
	{
		vector<ChimeraAlnType<_locr, _locl> > paired_seeds;
		pairChimeraSeeds(diagonalsegment_, paired_seeds,query_len);

		vector<DiagonalSeeds<_locr,_locl> > arm;
		for (size_t i = 0; i < paired_seeds.size(); i++)
		{
			arm.push_back(paired_seeds[i].arm1);
			if(paired_seeds[i].is_chimera)
			{
				if(paired_seeds[i].arm2.len > 0)
				{
					arm.push_back(paired_seeds[i].arm2);
				}
				
			}
		}
		if(arm.size()<6)
		{
			for (size_t i = 0; i < arm.size(); i++)
			{
				Hits<_locr,_locl> h = arm[i].hit_;
				local.push_back(local_match<_val> (h.seed_offset_, ref->data(h.subject_)));
				floating_sw(&query[h.seed_offset_],
						local.back(),
						padding[frame],
						ScoreMatrix::get().rawscore(VATParameters::gapped_xdrop),
						VATParameters::gap_open + VATParameters::gap_extend,
						VATParameters::gap_extend,
						transcript_buf,
						Traceback ());
				const int score = local.back().score_;
				std::pair<size_t,size_t> l = ReferenceSeqs<_val>::data_->local_position(h.subject_);
				matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
				anchored_transform(local.back(), l.second, h.seed_offset_);
				stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);
				to_source_space(local.back(), frame, dna_len);
				stat.inc(Statistics::SCORE_TOTAL, local.back().score_);
				stat.inc(Statistics::OUT_HITS);
			}

		}
		// cout<<"arm size = "<<arm.size()<<endl;


		// cout<<(int)h.query_<<"\t"<<(int)h.seed_offset_<<"\t"<<(int)h.subject_<<endl;
		// for (size_t i = 0; i < arm.size(); i++)
		// {
		// 	Hits<_locr,_locl> h = arm[i].hit_;
		// 	local.push_back(local_match<_val> (h.seed_offset_, ref->data(h.subject_)));
		// 	floating_sw(&query[h.seed_offset_],
		// 			local.back(),
		// 			padding[frame],
		// 			ScoreMatrix::get().rawscore(VATParameters::gapped_xdrop),
		// 			VATParameters::gap_open + VATParameters::gap_extend,
		// 			VATParameters::gap_extend,
		// 			transcript_buf,
		// 			Traceback ());
		// 	const int score = local.back().score_;
		// 	std::pair<size_t,size_t> l = ReferenceSeqs<_val>::data_->local_position(h.subject_);
		// 	matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
		// 	anchored_transform(local.back(), l.second, h.seed_offset_);
		// 	stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);
		// 	to_source_space(local.back(), frame, dna_len);
		// 	stat.inc(Statistics::SCORE_TOTAL, local.back().score_);
		// 	stat.inc(Statistics::OUT_HITS);
		// }

		/*
		for (size_t i = 0; i < paired_seeds.size(); i++)
		{
			Hits<_locr,_locl> h1 = paired_seeds[i].arm1.hit_;
			cout<<(int)h1.query_<<"\t"<<(int)h1.seed_offset_<<"\t"<<(int)h1.subject_<<endl;
			local.push_back(local_match<_val> (h1.seed_offset_, ref->data(h1.subject_)));
			floating_sw(&query[h1.seed_offset_],
					local.back(),
					padding[frame],
					ScoreMatrix::get().rawscore(VATParameters::gapped_xdrop),
					VATParameters::gap_open + VATParameters::gap_extend,
					VATParameters::gap_extend,
					transcript_buf,
					Traceback ());
			const int score = local.back().score_;
			std::pair<size_t,size_t> l = ReferenceSeqs<_val>::data_->local_position(h1.subject_);
			matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
			anchored_transform(local.back(), l.second, h1.seed_offset_);
			// stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);
			to_source_space(local.back(), frame, dna_len);
			// stat.inc(Statistics::SCORE_TOTAL, local.back().score_);
			// stat.inc(Statistics::OUT_HITS);

			if(paired_seeds[i].is_chimera)
			{
				Hits<_locr,_locl> h2 = paired_seeds[i].arm2.hit_;
				cout<<(int)h2.query_<<"\t"<<(int)h2.seed_offset_<<"\t"<<(int)h2.subject_<<endl;
				local.push_back(local_match<_val> (h2.seed_offset_, ref->data(h2.subject_)));
				floating_sw(&query[h2.seed_offset_],
						local.back(),
						padding[frame],
						ScoreMatrix::get().rawscore(VATParameters::gapped_xdrop),
						VATParameters::gap_open + VATParameters::gap_extend,
						VATParameters::gap_extend,
						transcript_buf,
						Traceback ());
				const int score = local.back().score_;
				std::pair<size_t,size_t> l = ReferenceSeqs<_val>::data_->local_position(h2.subject_);
				matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
				anchored_transform(local.back(), l.second, h2.seed_offset_);
				// stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);
				to_source_space(local.back(), frame, dna_len);
				// stat.inc(Statistics::SCORE_TOTAL, local.back().score_);
				// stat.inc(Statistics::OUT_HITS);
			}
				// stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);
				// to_source_space(local.back(), frame, dna_len);
				// stat.inc(Statistics::SCORE_TOTAL, local.back().score_);
				// stat.inc(Statistics::OUT_HITS);
			
		}*/
		

	}else
	{

		for (size_t i = 0; i < diagonalsegment_.size(); i++)
		{
			Hits<_locr,_locl> h = diagonalsegment_[i].hit_;
			local.push_back(local_match<_val> (h.seed_offset_, ref->data(h.subject_)));
			floating_sw(&query[h.seed_offset_],
					local.back(),
					padding[frame],
					ScoreMatrix::get().rawscore(VATParameters::gapped_xdrop),
					VATParameters::gap_open + VATParameters::gap_extend,
					VATParameters::gap_extend,
					transcript_buf,
					Traceback ());
			const int score = local.back().score_;
			std::pair<size_t,size_t> l = ReferenceSeqs<_val>::data_->local_position(h.subject_);
			matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
			anchored_transform(local.back(), l.second, h.seed_offset_);
			stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);
			to_source_space(local.back(), frame, dna_len);
			stat.inc(Statistics::SCORE_TOTAL, local.back().score_);
			stat.inc(Statistics::OUT_HITS);
		}
	}
	/*

	for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator i = begin; i != end; ++i) 
	{

		if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) {
			stat.inc(Statistics::DUPLICATES);
			continue;
		}
		local.push_back(local_match<_val> (i->seed_offset_, ref->data(i->subject_)));
		floating_sw(&query[i->seed_offset_],
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
	*/

	// cout<<"align_sequence 3"<<endl;

}

#endif /* ALIGN_SEQUENCE_H_ */
