

#ifndef ALIGN_SEQUENCE_H_
#define ALIGN_SEQUENCE_H_

#include <vector>
#include "../dp/floating_sw2.h"
#include "../dp/floating_sw.h"
#include "./pairChimera.h"
#include "./ungappedSeed.h"
#include "findwholeGenoSeeds.h"
#include "splicedunGappedSeed.h"
#include "dp_chain_chimera.h"
#include "chaining.h"
using std::vector;
template<typename _val>
bool cmp_seed_offset(const local_match<_val> &lhs, const local_match<_val> &rhs){
  return lhs.query_anchor_ < rhs.query_anchor_;
}

typedef score_vector2<int8_t> sv;

template<typename _val>
void align_group(vector<sequence<const _val> > &ref_suffix_set,
		 vector<sequence<const _val> > &ref_prefix_set,
		 const sequence<const _val> query,
		 int begin,
		 int end,
		 int band,
		 vector<local_match<_val> > &segments,
		 sv *traceback,
		 sv *Mv,
		 sv *Mh,
		 vector<char> &transcript_buf)
{
        int query_suffix_len, query_prefix_len, ref_suffix_len, ref_prefix_len, query_offset, i;
	const int query_len = query.length();
	query_offset = segments[begin].query_anchor_;
	query_suffix_len = query_len - query_offset;
	query_prefix_len = query_offset + 1;
	ref_suffix_len = query_suffix_len + band;
	ref_prefix_len = query_prefix_len + band;
	
	sequence<const _val> query_suffix (&query[query_offset], query_suffix_len, 0);
    sequence<const _val> query_prefix (&query[0], query_prefix_len, 0);
	ref_suffix_set.clear();
        ref_prefix_set.clear();
	
	for (i = begin; i <= end; i++){
	  sequence<const _val> ref_suffix (ReferenceSeqs<_val>::data_->window_infix_dp_right(segments[i].subject_position_, ref_suffix_len));
	  sequence<const _val> ref_prefix (ReferenceSeqs<_val>::data_->window_infix_dp_left(segments[i].subject_position_, ref_prefix_len));
	  ref_suffix_set.push_back(ref_suffix);
	  ref_prefix_set.push_back(ref_prefix);
	}
	
	floating_sw2(ref_suffix_set,
		     ref_prefix_set,
		     query_suffix,
		     query_prefix,
		     ScoreMatrix::get().rawscore(VATParameters::gapped_xdrop),
		     VATParameters::gap_open + VATParameters::gap_extend,
		     VATParameters::gap_extend,
		     begin,
		     band,
		     segments,
		     traceback,
		     Mv,
		     Mh,
		     transcript_buf,
		     int8_t());
}
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

	int GAP = -VATParameters::gap_open - VATParameters::gap_extend;
	unsigned band, max_slen, width;
	int k, local_temp_size, start;
	band = padding[frame];
	width = band * 2 + 2;

	vector<local_match<_val> > local_temp;
	local_temp.clear();
	
	if (VATParameters::aligner_mode == VATParameters::accuracy_model){
	  std::sort(begin, end, Hits<_locr,_locl>::cmp_normalized_subject);
	  
	  for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator i = begin; i != end; ++i) {
	    if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) {
	      stat.inc(Statistics::DUPLICATES);
	      continue;
	    }
	  
	    local_temp.push_back(local_match<_val> (i->seed_offset_, ref->data(i->subject_), i->subject_));
	  }
	}
	else {
	  for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator i = begin; i != end; ++i) {
	    local_temp.push_back(local_match<_val> (i->seed_offset_, ref->data(i->subject_), i->subject_));
          }
	}

	// sort for grouping
	std::sort(local_temp.begin(), local_temp.end(), cmp_seed_offset<_val>);
	
	// bottom local_match, which is useless
	local_temp.push_back(local_match<_val> (-1, ref->data(begin->subject_)));
	
	vector<sequence<const _val> > ref_suffix_set;
	vector<sequence<const _val> > ref_prefix_set;
	
	local_temp_size = local_temp.size() - 2;
	
	max_slen = max(query_len-local_temp[0].query_anchor_+band, local_temp[local_temp_size].query_anchor_+band+1);                                                                   
        max_slen += 2;

	sv *traceback = (sv*) malloc(max_slen * width * sizeof(sv));
        sv *Mv = (sv*)malloc(width * 3 * sizeof(sv));
        sv *Mh = (sv*)malloc(max_slen * 2 * sizeof(sv));
	
	unsigned a;
	// Mh
	memset(Mh, GAP, max_slen * 2 * sizeof(sv));
        Mh[0] = sv (0);   // mh
	// Dh
	for (a = 2; a <= band; a++){
          Mh[a*2] = sv (-1);
        }

        // traceback
	memset(traceback, BASE, max_slen * width * sizeof(sv));
	
	start = 0;
	for (k = 0; k <= local_temp_size; k++){
	  if (local_temp[k].query_anchor_ == local_temp[k+1].query_anchor_){
	    continue;
	  }
	  
	  // a group of hits range from local_temp[start] to local_temp[k]
	  align_group<_val>(ref_suffix_set, ref_prefix_set, query, start, k, band, local_temp, traceback, Mv, Mh, transcript_buf);
	  //count++;
	  start = k + 1;
	}
	//output_mutex.lock();
	//cout << "ha: " << count << endl;
	//output_mutex.unlock();

	for (k = 0; k <= local_temp_size; k++){
	  const int score = local_temp[k].score_;
          
	  local.push_back(local_temp[k]);
	  std::pair<size_t, size_t> l = ReferenceSeqs<_val>::data_->local_position(local.back().subject_position_);
	  matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
	  anchored_transform(local.back(), l.second, local.back().query_anchor_);
	  stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);
	  
	  //local.back().print(query, ref_seqs<_val>::get()[l.first], transcript_buf);
	  //cout << local.back().subject_begin_ << " " << local.back().subject_len_ << " " << local.back().query_begin_ << " " << local.back().query_len_ << endl;
	  
	  to_source_space(local.back(), frame, dna_len);
	  stat.inc(Statistics::SCORE_TOTAL, score);
	  stat.inc(Statistics::OUT_HITS);
	}

	
	free(traceback);
	free(Mv);
	free(Mh);

	//output_mutex.unlock();
	
	// vector<DiagonalSeeds<_locr,_locl> > diagonalsegment_;
	// // cout<<query<<endl;
	// for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator i = begin; i != end; ++i)
	// {
	// 	//(i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]
	// 	if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) 
	// 	{
	// 		stat.inc(Statistics::DUPLICATES);
	// 		continue;
	// 	}
	// 	std::pair<size_t, size_t> l = ref->local_position(i->subject_);
	// 	const _val* sbj = ref->data(i->subject_);
	// 	const _val* qry = &query[i->seed_offset_];
	// 	DiagonalSeeds<_locr,_locl> ds = ungappedSeeds<_val, _locr,_locl> (qry, sbj,(int)i->seed_offset_,(int)l.second,*i);
	// 	// if (ds.len >= VATParameters::gappex)
	// 	// {
	// 		diagonalsegment_.push_back(ds);
	// 	// }

		
	// }
	// // cout<<"seed num = "<<diagonalsegment_.size()<<endl;
	// if(VATParameters::chimera)
	// {
	// 	// vector<DiagonalSeeds<_locr,_locl> > chimera = diagonalsegment_;
	// 	int max_gap = 2000;
	// 	// cout<<"chimera = "<<endl;
	// 	vector<DiagonalSeeds<_locr,_locl> > chimera = findchimericSeeds(diagonalsegment_,max_gap);
	// 	// cout<<"chimera = "<<chimera.size()<<endl;
	// 	for (size_t i = 0; i < chimera.size(); i++)
	// 	{
	// 			Hits<_locr,_locl> h = chimera[i].hit_;
	// 			local.push_back(local_match<_val> (h.seed_offset_, ref->data(h.subject_)));
	// 			floating_sw(&query[h.seed_offset_],
	// 					local.back(),
	// 					padding[frame],
	// 					ScoreMatrix::get().rawscore(VATParameters::gapped_xdrop),
	// 					VATParameters::gap_open + VATParameters::gap_extend,
	// 					VATParameters::gap_extend,
	// 					transcript_buf,
	// 					Traceback ());
	// 			const int score = local.back().score_;
	// 			std::pair<size_t,size_t> l = ReferenceSeqs<_val>::data_->local_position(h.subject_);
	// 			matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
	// 			anchored_transform(local.back(), l.second, h.seed_offset_);
	// 			stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);
	// 			to_source_space(local.back(), frame, dna_len);
	// 			stat.inc(Statistics::SCORE_TOTAL, local.back().score_);
	// 			stat.inc(Statistics::OUT_HITS);

	// 	}

		
	// } else if(VATParameters::whole_genome)
	// {
	// 	/*
	// 	int maxDistance, int maxIndel, int maxLocalDistance
	// 	*/
		
	// 	vector<DiagonalSeeds<_locr,_locl> > whole_gen;
	// 	int max_gap = 50000; 
	// 	int maxIndel = 5; 
	// 	int maxLocalDistance = 40; 
	// 	whole_gen = ChainWGSSeeds(diagonalsegment_, max_gap,maxIndel,maxLocalDistance);
	// 	for (size_t i = 0; i < whole_gen.size(); i++)
	// 	{
	// 			Hits<_locr,_locl> h = whole_gen[i].hit_;
	// 			local.push_back(local_match<_val> (h.seed_offset_, ref->data(h.subject_)));
	// 			floating_sw(&query[h.seed_offset_],
	// 					local.back(),
	// 					padding[frame],
	// 					ScoreMatrix::get().rawscore(VATParameters::gapped_xdrop),
	// 					VATParameters::gap_open + VATParameters::gap_extend,
	// 					VATParameters::gap_extend,
	// 					transcript_buf,
	// 					Traceback ());
	// 			const int score = local.back().score_;
	// 			std::pair<size_t,size_t> l = ReferenceSeqs<_val>::data_->local_position(h.subject_);
	// 			matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
	// 			anchored_transform(local.back(), l.second, h.seed_offset_);
	// 			stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);
	// 			to_source_space(local.back(), frame, dna_len);
	// 			stat.inc(Statistics::SCORE_TOTAL, local.back().score_);
	// 			stat.inc(Statistics::OUT_HITS);

	// 	}
	// }
	// else if(VATParameters::spilce||VATParameters::circ)
	// {
	// 	// cout<<"spilce"<<endl;
	// 	vector<DiagonalSeeds<_locr,_locl> > spliced_seed;
	// 	const int max_gap = 250; 
	// 	spliced_seed = findSpliceSeeds(diagonalsegment_, max_gap);
	// 	for (size_t i = 0; i < spliced_seed.size(); i++)
	// 	{
	// 			Hits<_locr,_locl> h = spliced_seed[i].hit_;
	// 			local.push_back(local_match<_val> (h.seed_offset_, ref->data(h.subject_)));
	// 			floating_sw(&query[h.seed_offset_],
	// 					local.back(),
	// 					padding[frame],
	// 					ScoreMatrix::get().rawscore(VATParameters::gapped_xdrop),
	// 					VATParameters::gap_open + VATParameters::gap_extend,
	// 					VATParameters::gap_extend,
	// 					transcript_buf,
	// 					Traceback ());
	// 			const int score = local.back().score_;
	// 			std::pair<size_t,size_t> l = ReferenceSeqs<_val>::data_->local_position(h.subject_);
	// 			matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
	// 			anchored_transform(local.back(), l.second, h.seed_offset_);
	// 			stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);
	// 			to_source_space(local.back(), frame, dna_len);
	// 			stat.inc(Statistics::SCORE_TOTAL, local.back().score_);
	// 			stat.inc(Statistics::OUT_HITS);

	// 	}		
	// }
	// else
	// {
	// 	vector<DiagonalSeeds<_locr,_locl> > seeds;
	// 	int max_gap = 50000; 
	// 	int maxIndel = 0; 
	// 	int maxLocalDistance = 40; 
	// 	seeds = chainingSeeds(diagonalsegment_, max_gap,maxIndel,maxLocalDistance);

	// 	for (size_t i = 0; i < seeds.size(); i++)
	// 	{
	// 		Hits<_locr,_locl> h = seeds[i].hit_;
	// 		local.push_back(local_match<_val> (h.seed_offset_, ref->data(h.subject_)));
	// 		floating_sw(&query[h.seed_offset_],
	// 				local.back(),
	// 				padding[frame],
	// 				ScoreMatrix::get().rawscore(VATParameters::gapped_xdrop),
	// 				VATParameters::gap_open + VATParameters::gap_extend,
	// 				VATParameters::gap_extend,
	// 				transcript_buf,
	// 				Traceback ());
	// 		const int score = local.back().score_;
	// 		std::pair<size_t,size_t> l = ReferenceSeqs<_val>::data_->local_position(h.subject_);
	// 		matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
	// 		anchored_transform(local.back(), l.second, h.seed_offset_);
	// 		stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);
	// 		to_source_space(local.back(), frame, dna_len);
	// 		stat.inc(Statistics::SCORE_TOTAL, local.back().score_);
	// 		stat.inc(Statistics::OUT_HITS);
	// 	}
	// }
}

#endif /* ALIGN_SEQUENCE_H_ */
