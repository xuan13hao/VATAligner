

#ifndef ALIGN_UNGAPPED_H_
#define ALIGN_UNGAPPED_H_

// #include "../basic/DiagonalSegment.h"
#include "../basic/sequence.h"

/**
 * @brief xdropUngapped
 * 
 * @tparam _val 
 * @tparam _locr 
 * @tparam _locq 
 * @param query query
 * @param subject subject
 * @param seed_len seed len =  shape len
 * @param delta distance
 * @param len 
 * @return int score
 */

template<typename _val, typename _locr, typename _locq>
int xdropUngapped(const _val *query, const _val *subject, unsigned seed_len, unsigned &delta, unsigned &len)
{
	int score (0), st (0);
	unsigned n (0);
	delta = 0;
	// int count  = 0;
	const _val *q (query-1), *s (subject-1);
	const unsigned window_left = std::max(VATParameters::window, (unsigned)VATConsts::seed_anchor) - VATConsts::seed_anchor;
	while(score - st < VATParameters::xdrop
			&& delta < window_left
			&& *q != AlphabetSet<_val>::PADDING_CHAR
			&& *s != AlphabetSet<_val>::PADDING_CHAR)
	{
		// cout<<AlphabetAttributes<_val>::ALPHABET[*q];
		st += ScoreMatrix::get().letter_score(*q, mask_critical(*s));
		score = std::max(score, st);
		--q;
		--s;
		++delta;
	}

	q = query + seed_len;
	s = subject + seed_len;
	st = score;
	assert(seed_len >= VATConsts::seed_anchor);
	const unsigned window_right = std::max(VATParameters::window, seed_len - VATConsts::seed_anchor) - (seed_len - VATConsts::seed_anchor);
	while(score - st < VATParameters::xdrop
			&& n < window_right
			&& *q != AlphabetSet<_val>::PADDING_CHAR
			&& *s != AlphabetSet<_val>::PADDING_CHAR)
	{

		st += ScoreMatrix::get().letter_score(*q, mask_critical(*s));
		score = std::max(score, st);
		++q;
		++s;
		++n;
	}

	for(unsigned i=0;i<seed_len;++i)
		score += ScoreMatrix::get().letter_score(query[i], mask_critical(subject[i]));

	len = delta + n + seed_len;
	// return DiagonalSegment(qa - delta, sa - delta, len + delta, score)
	return score;
}
template<typename _val>
DiagonalSeeds xdrop_ungapped(const _val *query, const _val *subject, int qa, int sa, hit& h)
{
	unsigned delta,len;
	const int xdrop = 23;
	int score(0), st(0), n=1;
	delta = 0;
	const _val *q(query - 1), *s(subject - 1);
	// cout<<"qa = "<<qa<<", sa = "<<qa<<endl;
	while (score - st < xdrop
		&& *q != AlphabetSet<_val>::PADDING_CHAR
		&& *s != AlphabetSet<_val>::PADDING_CHAR)
	{
		// cout<<"ql = "<<AlphabetAttributes<_val>::ALPHABET[*q]<<",sl = "<<AlphabetAttributes<_val>::ALPHABET[*s]<<endl;
		st += ScoreMatrix::get().letter_score(*q, mask_critical(*s));
		if (st > score) {
			score = st;
			delta = n;
		}
		--q;
		--s;
		++n;
	}

	q = query;
	s = subject;
	st = score;
	n = 1;
	len = 0;
	while (score - st < xdrop
		&& *q != AlphabetSet<_val>::PADDING_CHAR
		&& *s != AlphabetSet<_val>::PADDING_CHAR)
	{
		// cout<<"qr = "<<AlphabetAttributes<_val>::ALPHABET[*q]<<",sr = "<<AlphabetAttributes<_val>::ALPHABET[*s]<<endl;
		st += ScoreMatrix::get().letter_score(*q, mask_critical(*s));
		if (st > score) {
			score = st;
			len = n;
		}
		++q;
		++s;
		++n;
	}
	// cout<<"qa = "<<qa<<", sa = "<<qa<<endl;
	len += delta;
	//int query_pos, int subject_pos, int len, int score
	// cout<<"i = "<<qa - delta<<",j ="<<sa - delta<<", len = "<<len + delta<<", score = "<<score<<", delta = "<<delta<<", len = "<<len<<endl;
	return DiagonalSeeds(qa - delta, sa - delta, len, score,h);
}
/*
template<typename _val, typename _locr, typename _locq>
DiagonalSegment xdrop_ungapped(const _val *query, const _val *subject, int qa, int sa)
{
	const int xdrop = VATParameters::xdrop;
	int score = 0, st = 0;
	int n = 1, delta = 0, len = 0;

	int q = qa - 1, s = sa - 1;
	const _val *ql (query-1), *sl (subject-1);
	while (score - st < xdrop
			&& *ql != AlphabetSet<_val>::PADDING_CHAR
			&& *sl != AlphabetSet<_val>::PADDING_CHAR)
	{
		st += ScoreMatrix::get().letter_score(*ql, mask_critical(*sl));
		if (st > score) 
		{
			score = st;
			delta = n;
		}
		--q;
		--s;
		--ql;
		--sl;
		++n;
	}

	q = qa;
	s = sa;
	st = score;
	n = 1;
	while (score - st < xdrop
			&& *ql != AlphabetSet<_val>::PADDING_CHAR
			&& *sl != AlphabetSet<_val>::PADDING_CHAR)
	{
		st += ScoreMatrix::get().letter_score(*ql, mask_critical(*sl));
		if (st > score) {
			score = st;
			len = n;
		}
		++q;
		++s;
		++ql;
		++sl;
		++n;
	}
	return DiagonalSegment(qa - delta, sa - delta, len + delta, score);
}
*/
#endif /* ALIGN_UNGAPPED_H_ */
