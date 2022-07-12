

#ifndef ALIGN_UNGAPPED_H_
#define ALIGN_UNGAPPED_H_
/**
 * 
 * GT.....AG intron
 */

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
	int count  = 0;
	const _val *q (query-1), *s (subject-1);
	const unsigned window_left = std::max(VATParameters::window, (unsigned)VATConsts::seed_anchor) - VATConsts::seed_anchor;
	while(score - st < VATParameters::xdrop
			&& delta < window_left
			&& *q != AlphabetSet<_val>::PADDING_CHAR
			&& *s != AlphabetSet<_val>::PADDING_CHAR)
	{
		
		if (VATParameters::is_spilced)
		{
			int spliced_reward = 0;
			// cout<<"query = "<<AlphabetAttributes<_val>::ALPHABET[*--q] <<",subject = "<<AlphabetAttributes<_val>::ALPHABET[*s]<<endl;
			if ((AlphabetAttributes<_val>::ALPHABET[*q] == 'T' && AlphabetAttributes<_val>::ALPHABET[*s] == 'T')
			&& (AlphabetAttributes<_val>::ALPHABET[*--q] == 'G' && AlphabetAttributes<_val>::ALPHABET[*--s] == 'G')
			&& (*--q != AlphabetSet<_val>::PADDING_CHAR)&&(*--s != AlphabetSet<_val>::PADDING_CHAR)
			)
			{
				// cout<<"count = "<<(count++)<<endl;
				spliced_reward = 12;
			}
			
			int letter_score = ScoreMatrix::get().letter_score(*q, mask_critical(*s)) + spliced_reward;
			st += letter_score;
			score = std::max(score, st);
			--q;
			--s;
			++delta;
		}
		
		
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
		if (VATParameters::is_spilced)
		{
			int spliced_reward = 0;
			// cout<<"query = "<<AlphabetAttributes<_val>::ALPHABET[*--q] <<",subject = "<<AlphabetAttributes<_val>::ALPHABET[*s]<<endl;
			if ((AlphabetAttributes<_val>::ALPHABET[*q] == 'T'&& AlphabetAttributes<_val>::ALPHABET[*s] == 'T')
			&& (AlphabetAttributes<_val>::ALPHABET[*++q] == 'G' && AlphabetAttributes<_val>::ALPHABET[*++s] == 'G')
			&& (*++q != AlphabetSet<_val>::PADDING_CHAR)&& (*++s != AlphabetSet<_val>::PADDING_CHAR)
			)
			{
				// cout<<"count = "<<(count++)<<endl;
				spliced_reward = 12;
			}
			
			int letter_score = ScoreMatrix::get().letter_score(*q, mask_critical(*s)) + spliced_reward;
			st += letter_score;
			score = std::max(score, st);
			--q;
			--s;
			++delta;
		}
		


		st += ScoreMatrix::get().letter_score(*q, mask_critical(*s));
		score = std::max(score, st);
		++q;
		++s;
		++n;
	}

	for(unsigned i=0;i<seed_len;++i)
		score += ScoreMatrix::get().letter_score(query[i], mask_critical(subject[i]));

	len = delta + n + seed_len;
	return score;
}

#endif /* ALIGN_UNGAPPED_H_ */
