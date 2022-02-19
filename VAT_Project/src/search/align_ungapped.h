

#ifndef ALIGN_UNGAPPED_H_
#define ALIGN_UNGAPPED_H_

template<typename _val, typename _locr, typename _locq>
int xdrop_ungapped(const _val *query, const _val *subject, unsigned seed_len, unsigned &delta, unsigned &len)
{
	int score (0), st (0);
	unsigned n (0);
	delta = 0;

	const _val *q (query-1), *s (subject-1);
	const unsigned window_left = std::max(VATParameters::window, (unsigned)VATConsts::seed_anchor) - VATConsts::seed_anchor;
	while(score - st < VATParameters::xdrop
			&& delta < window_left
			&& *q != String_set<_val>::PADDING_CHAR
			&& *s != String_set<_val>::PADDING_CHAR)
	{
		st += score_matrix::get().letter_score(*q, mask_critical(*s));
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
			&& *q != String_set<_val>::PADDING_CHAR
			&& *s != String_set<_val>::PADDING_CHAR)
	{
		st += score_matrix::get().letter_score(*q, mask_critical(*s));
		score = std::max(score, st);
		++q;
		++s;
		++n;
	}

	for(unsigned i=0;i<seed_len;++i)
		score += score_matrix::get().letter_score(query[i], mask_critical(subject[i]));

	len = delta + n + seed_len;
	return score;
}

#endif /* ALIGN_UNGAPPED_H_ */
