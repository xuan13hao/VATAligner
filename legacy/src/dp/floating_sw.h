

#ifndef FLOATING_SW_H_
#define FLOATING_SW_H_

#include "../basic/Hits.h"
#include "scalar_dp_matrix.h"
#include "../tools/direction.h"
#include "scalar_traceback.h"

template<typename _val, typename _dir, typename _score, typename _traceback>
local_match<_val> floating_sw_dir(const _val *query, const _val *subject, int band, _score xdrop, _score gap_open, _score gap_extend, vector<char> &transcript_buf)
{
	using std::max;

	_score max_score = 0, column_max = 0;
	int j = 0, i_max = -1, j_best = -1, i_best = -1;
	Scalar_dp_matrix<_score,_traceback> mtx (band);
	const _val *x = query, *y = subject;
	
	while(*y != AlphabetSet<_val>::PADDING_CHAR && max_score - column_max < xdrop) 
	{
		typename Scalar_dp_matrix<_score,_traceback>::Column_iterator it = mtx.column(j, i_max);
		if(get_dir(x, it.row(), _dir()) == AlphabetSet<_val>::PADDING_CHAR)
			break;
		_score vgap = Scalar_dp_matrix<_score,_traceback>::NEG_MIN;
		if(get_dir(x, i_max+1, _dir()) == AlphabetSet<_val>::PADDING_CHAR) {
			column_max = std::numeric_limits<_score>::min();
		} else {
			++i_max;
			column_max += ScoreMatrix::get().letter_score(mask_critical(*y), get_dir(x, i_max, _dir()));
		}
		// cout<<(int)get_dir(x, it.row(), _dir())<<endl;
		for(; it.valid() && get_dir(x, it.row(), _dir()) != AlphabetSet<_val>::PADDING_CHAR; ++it) 
		{
			const _score match_score = ScoreMatrix::get().letter_score(mask_critical(*y), get_dir(x, it.row(), _dir()));
			// cout<<"y= "<<(int)mask_critical(*y)<<", x = "<<(int)get_dir(x, it.row(), _dir())<<endl;
			const _score s = max(max(it.diag() + match_score, vgap), it.hgap_in());
			if(s > column_max) 
			{
				column_max = s;
				i_max = it.row();
			}
			const _score open = s - gap_open;
			vgap = max(vgap - gap_extend, open);
			it.hgap_out() = max(it.hgap_in() - gap_extend, open);
			it.score() = s;
		}

		if(column_max > max_score) {
			max_score = column_max;
			j_best = j;
			i_best = i_max;
		}
		y = inc_dir(y, _dir());
		++j;
	}

	return traceback<_val,_dir,_score>(query, subject, mtx.score_buffer(), band, gap_open, gap_extend, j_best, i_best, max_score, transcript_buf);
}

template<typename _val, typename _score, typename _traceback>
void floating_sw(const _val *query, local_match<_val> &segment, int band, _score xdrop, _score gap_open, _score gap_extend, vector<char> &transcript_buf, const _traceback& = Score_only (), const _score& = int())
{
	// cout<<"floating_sw"<<endl;
	segment += floating_sw_dir<_val,Right,_score,_traceback>(query, segment.subject_, band, xdrop, gap_open, gap_extend, transcript_buf);
	const local_match<_val> left (floating_sw_dir<_val,Left,_score,_traceback>(query, segment.subject_, band, xdrop, gap_open, gap_extend, transcript_buf));
	if(left.query_len_ > 0) {
		segment -= left;
		segment.query_begin_--;
		segment.subject_begin_--;
		const _val q = *query, s = mask_critical(*segment.subject_);
		segment.score_ -= ScoreMatrix::get().letter_score(q, s);
		if(q == s)
			segment.identities_--;
		else
			segment.mismatches_--;
		segment.len_--;
		segment.subject_len_--;
		segment.query_len_--;
	}
}

#endif /* FLOATING_SW_H_ */
