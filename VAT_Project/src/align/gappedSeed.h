#ifndef __GAPPEDSEED_H__
#define __GAPPEDSEED_H__

template<typename _val>
vector<DiagonalSeeds> getSeeds(vector<Segment<_val> > &matches, int query_len, int query)
{
	typename vector<Segment<_val> >::iterator it = matches.begin();
	const int min_raw_score = ScoreMatrix::get().rawscore(VATParameters::min_bit_score == 0
			? ScoreMatrix::get().bitscore(VATParameters::max_evalue, ref_header.letters, query_len) : VATParameters::min_bit_score);	
	vector<DiagonalSeeds> rawseeds;
	while(it < matches.end()) 
	{
		int q_id,s_id,q_pos,s_pos,score,length;
		q_id = query;
		s_id = it->subject_id_;
		length = it->traceback_->len_;
		score =  it->score_;
		q_pos = it->traceback_->query_begin_;
		s_pos = it->traceback_->subject_begin_;
		DiagonalSeeds ds(q_pos,s_pos,length,score,q_id,s_id);
		rawseeds.push_back(ds);
		++it;
	}
	return rawseeds;
}

#endif // __GAPPEDSEED_H__