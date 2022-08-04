#ifndef __UNGAPPEDSEEDS_H__
#define __UNGAPPEDSEEDS_H__


#include "../basic/DiagonalSegment.h"
/**
 * 
 * AlphabetSet<_val>::PADDING_CHAR is a mark of the end of a sequence 
 * AlphabetAttributes<_val>::ALPHABET[s.data_[i]];
 * 
 * GT.....AG intron
 * 
 * AlphabetAttributes<_val>::ALPHABET[x] = T;
 * AlphabetAttributes<_val>::ALPHABET[x-1] = G;
 * 
 * AlphabetAttributes<_val>::ALPHABET[x] = A;
 * AlphabetAttributes<_val>::ALPHABET[x+1] = G;
 */

template<typename _val, typename _locr, typename _locq>
DiagonalSeeds xdrop_ungapped(const _val &query, const _val &subject, int qa, int sa)
{
	const int xdrop = VATParameters::xdrop;
	int score = 0, st = 0;
	int n = 1, delta = 0, len = 0;

	int q = qa - 1, s = sa - 1;
	_val ql, sl;
	while (score - st < xdrop
		&& (ql = query[q]) != AlphabetSet<_val>::PADDING_CHAR
		&& (sl = subject[s]) != AlphabetSet<_val>::PADDING_CHAR)
	{
		st += ScoreMatrix::get().letter_score(ql, mask_critical(sl));
		// st += score_matrix(ql, sl);
		if (st > score) {
			score = st;
			delta = n;
		}
		--q;
		--s;
		++n;
	}

	q = qa;
	s = sa;
	st = score;
	n = 1;
	while (score - st < xdrop
		&& (ql = query[q]) != AlphabetSet<_val>::PADDING_CHAR
		&& (sl = subject[s]) != AlphabetSet<_val>::PADDING_CHAR)
	{
		st += ScoreMatrix::get().letter_score(ql, mask_critical(sl));
		if (st > score) {
			score = st;
			len = n;
		}
		++q;
		++s;
		++n;
	}
	return DiagonalSeeds(qa - delta, sa - delta, len + delta, score);
}
/*

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
















*/

#endif // __UNGAPPEDSEEDS_H__