

#ifndef ALIGN_SEQUENCE_H_
#define ALIGN_SEQUENCE_H_

#include <vector>
#include <chrono>
#include "../dp/floating_sw.h"
#include "./pairChimera.h"
#include "./ungappedSeed.h"
#include "findwholeGenoSeeds.h"
#include "splicedunGappedSeed.h"
#include "dp_chain_chimera.h"
#include "chaining.h"
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
	// std::cout << "align_sequence" <<std::endl;
	std::sort(begin, end, Hits<_locr,_locl>::cmp_normalized_subject);
	const unsigned q_num (begin->query_);
	const sequence<const _val> query (QuerySeqs<_val>::get()[q_num]);
	const unsigned frame = q_num % query_contexts();
	const unsigned query_len = query.length();
	padding[frame] = VATParameters::read_padding<_val>(query_len);
	const SequenceSet<_val> *ref = ReferenceSeqs<_val>::data_;

	if(VATParameters::chimera)
	{
		vector<DiagonalSeeds<_locr,_locl> > diagonalsegment_;
		for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator i = begin; i != end; ++i)
		{
			//(i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]
			if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) 
			{
				stat.inc(Statistics::DUPLICATES);
				continue;
			}
			std::pair<size_t, size_t> l = ref->local_position(i->subject_);
			const _val* sbj = ref->data(i->subject_);
			const _val* qry = &query[i->seed_offset_];
			DiagonalSeeds<_locr,_locl> ds = ungappedSeeds<_val, _locr,_locl> (qry, sbj,(int)i->seed_offset_,(int)l.second,*i);
			if (ds.len >= 10)
			{
				diagonalsegment_.push_back(ds);
			}

			
		}
		// vector<DiagonalSeeds<_locr,_locl> > chimera = diagonalsegment_;
		int max_gap = 2000;
		// cout<<"chimera = "<<endl;
		vector<DiagonalSeeds<_locr,_locl> > chimera = findchimericSeeds(diagonalsegment_,max_gap);
		// cout<<"chimera = "<<chimera.size()<<endl;
		for (size_t i = 0; i < chimera.size(); i++)
		{
				Hits<_locr,_locl> h = chimera[i].hit_;
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

		
	} else if(VATParameters::whole_genome)
	{

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
			to_source_space(local.back(), frame, dna_len);
			stat.inc(Statistics::SCORE_TOTAL, local.back().score_);
			stat.inc(Statistics::OUT_HITS);
		}
	}
else if(VATParameters::dna_homology)
	{

		for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator i = begin; i != end; ++i) 
		{
			if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) 
			{
				stat.inc(Statistics::DUPLICATES);
				continue;
			}
			local.push_back(local_match<_val> (i->seed_offset_, ref->data(i->subject_)));
			floating_sw(&query[i->seed_offset_],
					local.back(),
					padding[frame],
					VATParameters::gapped_xdrop,
					VATParameters::gap_open,
					VATParameters::gap_extend,
					transcript_buf,
					Traceback ());
			const int score = local.back().score_;
			std::pair<size_t,size_t> l = ReferenceSeqs<_val>::data_->local_position(i->subject_);
			matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
			anchored_transform(local.back(), l.second, i->seed_offset_);
			stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);
			to_source_space(local.back(), frame, dna_len);
			stat.inc(Statistics::SCORE_TOTAL, local.back().score_);
			stat.inc(Statistics::OUT_HITS);
		}
	}
	else if(VATParameters::spilce||VATParameters::circ)
	{
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
			if (ds.len >= 10)
			{
				diagonalsegment_.push_back(ds);
			}
		}
		vector<DiagonalSeeds<_locr,_locl> > spliced_seed;
		const int max_gap = 250; 
		spliced_seed = findSpliceSeeds(diagonalsegment_, max_gap);
		for (size_t i = 0; i < spliced_seed.size(); i++)
		{
				Hits<_locr,_locl> h = spliced_seed[i].hit_;
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
	else
	{
		vector<DiagonalSeeds<_locr,_locl> > diagonalsegment_;
		for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator i = begin; i != end; ++i)
		{
			// (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]
			if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) 
			{
				stat.inc(Statistics::DUPLICATES);
				continue;
			}
			std::pair<size_t, size_t> l = ref->local_position(i->subject_);
			const _val* sbj = ref->data(i->subject_);
			const _val* qry = &query[i->seed_offset_];
			DiagonalSeeds<_locr,_locl> ds = ungappedSeeds<_val, _locr,_locl> (qry, sbj,(int)i->seed_offset_,(int)l.second,*i);
			if (ds.len >= 5)
			{
				diagonalsegment_.push_back(ds);
			}

			
		}
		vector<DiagonalSeeds<_locr,_locl> > seeds;
		int max_gap = 50000; 
		int maxIndel = 0; 
		int maxLocalDistance = 40;
		 
		seeds = chainingSeeds(diagonalsegment_, max_gap,maxIndel,maxLocalDistance);

		for (size_t i = 0; i < seeds.size(); i++)
		{
			Hits<_locr,_locl> h = seeds[i].hit_;
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
}

#endif /* ALIGN_SEQUENCE_H_ */
