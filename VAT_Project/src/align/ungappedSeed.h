#ifndef __UNGAPPEDSEED_H__
#define __UNGAPPEDSEED_H__
#include <algorithm>
#include "../basic/Seed.h"
using namespace std;
template<typename _val, typename _locr, typename _locq>
DiagonalSeeds<_locr, _locq> ungappedSeeds(const _val *query, const _val *subject, int qa, int sa, Hits<_locr, _locq>& h)
{
	int sbj_, qry_;
	qry_=h.query_;
	string q_;
	string r_;
	std::pair<size_t,size_t> l = ReferenceSeqs<_val>::data_->local_position(h.subject_);
	sbj_ = l.first;
	unsigned delta,len;
	int score(0), st(0), n=1;
	int match(0);
	int tmp_len(0);
	int mismatch(0);
	delta = 0;
	const _val *q(query - 1), *s(subject - 1);
	while (score - st < VATParameters::gapped_xdrop
		&& *q != AlphabetSet<_val>::PADDING_CHAR
		&& *s != AlphabetSet<_val>::PADDING_CHAR)
	{
		q_.push_back(AlphabetAttributes<_val>::ALPHABET[*q]);
		r_.push_back(AlphabetAttributes<_val>::ALPHABET[*s]);
		// st += ScoreMatrix::get().letter_score(*q, mask_critical(*s));
		// if (*q == mask_critical(*s)) {
		// 	tmp_len = 1;
		// 	match++;
		// } else {
		// 	tmp_len = 0;
		// 	mismatch++;
		// }
		// st +=tmp_len;
		if (st > score) {
			score = st;
			delta = n;
		}
		--q;
		--s;
		++n;
	}
	reverse(q_.begin(), q_.end());
	reverse(r_.begin(), r_.end());
	q = query;
	s = subject;
	st = score;
	n = 1;
	len = 0;
	while (score - st < VATParameters::gapped_xdrop
		&& *q != AlphabetSet<_val>::PADDING_CHAR
		&& *s != AlphabetSet<_val>::PADDING_CHAR)
	{
		q_.push_back(AlphabetAttributes<_val>::ALPHABET[*q]);
		r_.push_back(AlphabetAttributes<_val>::ALPHABET[*s]);
		st += ScoreMatrix::get().letter_score(*q, mask_critical(*s));
		// if (*q == mask_critical(*s)) {
		// 	tmp_len = 1;
		// 	match++;
		// } else {
		// 	tmp_len = 0;
		// 	mismatch++;
		// }
		// st +=tmp_len;
		if (st > score) {
			score = st;
			len = n;
		}
		++q;
		++s;
		++n;
	}
	// cout<<"q = "<<q_<<", r = "<<r_<<endl;
	len += delta;

	// if (len > VATParameters::seed_len) {
    //     len = VATParameters::seed_len;
    //     q_.resize(len);
    //     r_.resize(len);
    // }
	//int query_pos, int subject_pos, int len, int score
	// cout<<"i = "<<qa - delta<<",j ="<<sa - delta<<", len = "<<len + delta<<", score = "<<score<<", delta = "<<delta<<", len = "<<len<<endl;
	return DiagonalSeeds<_locr, _locq>(qa - delta, sa - delta, len, score,match,mismatch,h,qry_,sbj_,q_,r_);
}

#endif // __UNGAPPEDSEED_H__