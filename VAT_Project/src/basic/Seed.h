#ifndef __SEED_H__
#define __SEED_H__

#include <string> // for string class
#include <iostream>
#include <algorithm>
#include "Hits.h"

using namespace std;
/**
 * DiagonalSeeds i, j, len, score
 * int query_pos, int subject_pos, int len, int score, string q ,s
*/
template<typename _locr, typename _locl>
class DiagonalSeeds
{
    public:
	DiagonalSeeds() :
		len(0)
	{}
	DiagonalSeeds(int query_pos, int subject_pos, int len, int score) :
		i(query_pos),
		j(subject_pos),
		len(len),
		score(score)
	{
		
	}
	DiagonalSeeds(int query_pos, int subject_pos, int len, int score, int match,int mismatch,Hits<_locr,_locl>& h, int& q_, int& s_,string& q, string& s) :
	i(query_pos),
	j(subject_pos),
	len(len),
	score(score),
	hit_(h),
	qry_id(s_),
	sbj_id(q_),
	qry_str(q),
	match(match),
	mismatch(mismatch),
	sbj_str(s)
	{
		
	}

	bool empty() const
	{
		return len == 0;
	}
	int subject_last() const
	{
		return j + len - 1;
	}
	int query_last() const
	{
		return i + len - 1;
	}
	int subject_end() const
	{
		return j + len;
	}
	int query_end() const
	{
		return i + len;
	}
	int diag() const
	{
		return i - j;
	}
	DiagonalSeeds transpose() const
	{
		return DiagonalSeeds(j, i, len, score);
	}
	int partial_score(int diff) const
	{
		return score*std::max(len - diff, 0) / len;
	}
	bool operator<=(const DiagonalSeeds &rhs) const
	{
		return i + len <= rhs.i && j + len <= rhs.j;
	}

	DiagonalSeeds& operator=(const DiagonalSeeds &rhs)
	{
		
		i = rhs.i;
		j = rhs.j;
		hit_ = rhs.hit_;
		qry_id = rhs.qry_id;
		sbj_id = rhs.sbj_id;
		score = rhs.score;
		len = rhs.len;
		qry_str = rhs.qry_str;
		sbj_str = rhs.sbj_str;
		match = rhs.match;
		mismatch = rhs.mismatch;
		return *this;
	}
	bool operator==(const DiagonalSeeds &rhs) const
	{
		return i == rhs.i && j == rhs.j && len == rhs.len;
	}
	static bool cmp_subject(const DiagonalSeeds &x, const DiagonalSeeds &y)
	{
		return x.j < y.j || (x.j == y.j && x.i < y.i);
	}
	static bool cmp_score(const DiagonalSeeds& x, const DiagonalSeeds& y)
	{
		return x.score > y.score;
	}
	static bool cmp_query_begin(const DiagonalSeeds &x, const DiagonalSeeds &y)
	{
		return x.i < y.i;
	}
	static bool cmp_subject_end(const DiagonalSeeds &x, const DiagonalSeeds &y)
	{
		return x.subject_end() < y.subject_end();
	}

	bool isChimericMapping(const DiagonalSeeds& other) const {
        // Check if seeds correspond to different reference sequences
        bool differentReference = (i != other.i);
		int MAX_DISTANCE_IN_REFERENCE = 5;
        // Check if seeds are too far away in a single reference sequence
        bool farAwayInReference = (std::abs(static_cast<int>(j - other.j)) > MAX_DISTANCE_IN_REFERENCE);

        // Determine if seeds induce chimeric mapping
        return differentReference || farAwayInReference;
    }
	friend std::ostream& operator<<(std::ostream &s, const DiagonalSeeds &d)
	{
		s << "i=" << d.i << " j=" << d.j << " l=" << d.len << " score=" << d.score<<" match = "<<d.match<<" mismatch = "<<d.mismatch;
		return s;
	}
	int i, j, len, score;//query_pos, subject_pos
	string qry_str;
	string sbj_str;
	int qry_id, sbj_id;
	int match,mismatch;
	Hits<_locr,_locl> hit_;
};


#endif // __SEED_H__