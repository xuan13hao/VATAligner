

#ifndef SCALAR_TRACEBACK_H_
#define SCALAR_TRACEBACK_H_

template<typename _score>
struct Scalar_traceback_matrix
{
	Scalar_traceback_matrix(const Growing_buffer<_score> &data, int band):
		data_ (data),
		band_ (band)
	{ }
	int operator()(int col, int row) const
	{ return data_.column(col+1)[row - (data_.center(col+1)-band_)]; }
	bool in_band(int col, int row) const
	{ return row >= data_.center(col+1)-band_ && row <= data_.center(col+1)+band_ && row >= 0 && col >= 0; }
	void print(int col, int row) const
	{
		for(unsigned j=0;j<=row;++j) {
			for(unsigned i=0;i<=col;++i)
				printf("%4i", in_band(i, j) ? this->operator()(i, j) : 0);
			printf("\n");
		}
	}
private:
	const Growing_buffer<_score> &data_;
	const int band_;
};

template<typename _score>
bool have_vgap(const Scalar_traceback_matrix<_score> &dp,
		int i,
		int j,
		int gap_open,
		int gap_extend,
		int &l)
{
	int score = dp(i, j);
	l = 1;
	--j;
	while(dp.in_band(i, j)) {
		if(score == dp(i, j) - gap_open - (l-1)*gap_extend)
			return true;
		--j;
		++l;
	}
	return false;
}

template<typename _score>
bool have_hgap(const Scalar_traceback_matrix<_score> &dp,
		int i,
		int j,
		int gap_open,
		int gap_extend,
		int &l)
{
	int score = dp(i, j);
	l = 1;
	--i;
	while(dp.in_band(i, j)) {
		if(score == dp(i, j) - gap_open - (l-1)*gap_extend)
			return true;
		--i;
		++l;
	}
	return false;
}

template<typename _val, typename _dir, typename _score>
local_match<_val> traceback(const _val *query,
		const _val *subject,
		const Growing_buffer<_score> &scores,
		int band,
		int gap_open,
		int gap_extend,
		int i,
		int j,
		int score,
		vector<char> &transcript_buf)
{
	if(i == -1)
		return local_match<_val> (0);
	Scalar_traceback_matrix<_score> dp (scores, band);
	//dp.print(i, j);

	local_match<_val> l;
	l.query_len_ = j + 1;
	l.subject_len_ = i + 1;
	l.query_begin_ = 0;
	l.subject_begin_ = 0;
	l.score_ = score;
	Edit_transcript transcript (transcript_buf);

	int gap_len;

	while(i>0 || j>0) 
	{
		const _val lq = get_dir(query, j, _dir()), ls = mask_critical(get_dir(subject, i, _dir()));
		const int match_score = ScoreMatrix::get().letter_score(lq, ls);
		// printf("i=%i j=%i score=%i subject=%c query=%c\n",i,j,dp(i, j),AlphabetAttributes<_val>::ALPHABET[ls],AlphabetAttributes<_val>::ALPHABET[lq]);

		if(dp(i, j) == match_score + dp(i-1, j-1)) 
		{
			if(lq == ls)
				++l.identities_;
			else
				++l.mismatches_;
			--i;
			--j;
			++l.len_;
			transcript_buf.push_back(op_match);
		} else if (have_hgap(dp, i, j, gap_open, gap_extend, gap_len)) {
			++l.gap_openings_;
			l.len_ += gap_len;
			i -= gap_len;
			transcript_buf.insert(transcript_buf.end(), gap_len ,op_deletion);
		} else if (have_vgap(dp, i, j, gap_open, gap_extend, gap_len)) {
			++l.gap_openings_;
			l.len_ += gap_len;
			j -= gap_len;
			transcript_buf.insert(transcript_buf.end(), gap_len ,op_insertion);
		} else 
		{
			// dp.print(i,j);
			throw std::runtime_error("Traceback error.");
		}
	}

	if(get_dir(query, 0, _dir()) == mask_critical(get_dir(subject, 0, _dir())))
		++l.identities_;
	else
		++l.mismatches_;
	++l.len_;
	transcript_buf.push_back(op_match);
	//printf("len=%i\n",l.len_);
	l.transcript_right_ = transcript.set_end(transcript_buf);
	return l;
}

template<typename _val, typename _dir, typename _score>
local_match<_val> traceback(const _val *query,
		const _val *subject,
		const Double_buffer<_score> &scores,
		int band,
		int gap_open,
		int gap_extend,
		int i,
		int j,
		int score)
{ return local_match<_val> (score); }

#endif /* SCALAR_TRACEBACK_H_ */
