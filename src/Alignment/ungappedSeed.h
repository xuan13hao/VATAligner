#ifndef __UNGAPPEDSEED_H__
#define __UNGAPPEDSEED_H__
#include <algorithm>
#include "../Commons/Seed.h"
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
	//DiagonalSeeds(int query_pos, int subject_pos, int len, int score,Hits<_locr,_locl>& h, bool splice, bool back_splice)
    bool i_spjunction = false;
    if (r_.substr(0, 2) == "GT" && r_.substr(r_.length() - 2, 2) == "AG") 
	{
            i_spjunction=  true;
	}
	bool BackSplicedJunction = false;
    if (len < 8) {
        BackSplicedJunction = false; // Seed is too short to be a candidate junction
    }
    // Check for splice site motifs
    size_t donorPos = q_.find("GT");
    size_t acceptorPos = q_.find("AG");

    // Check if both donor and acceptor motifs are present
    if (donorPos == std::string::npos || acceptorPos == std::string::npos) {
        BackSplicedJunction = false; // Seed does not contain both donor and acceptor motifs
    }

    // Check orientation
    if (donorPos > acceptorPos) {
        BackSplicedJunction = true; // Seed has the reverse orientation of splice sites
    }

	return DiagonalSeeds<_locr, _locq>(qa - delta, sa - delta, len, score,h,i_spjunction,BackSplicedJunction);
	// return DiagonalSeeds<_locr, _locq>(qa - delta, sa - delta, len, score,match,mismatch,h,qry_,sbj_,q_,r_);
}

#endif // __UNGAPPEDSEED_H__