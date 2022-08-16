#ifndef __DIAGONALSEGMENT_H__
#define __DIAGONALSEGMENT_H__

#include <string> // for string class
#include <iostream>
#include <algorithm>
#include "interval.h"
#include "Hits.h"
using namespace std;
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

	DiagonalSeeds(int query_pos, int subject_pos, int len, int score, hit& h, string& q_, string& s_) :
		i(query_pos),
		j(subject_pos),
		len(len),
		score(score),
		hit_(h),
		sbj_(s_),
		qry_(q_)
	{
		
	}
	bool empty() const
	{
		return len == 0;
	}
	interval query_range() const
	{
		return interval(i, i + len);
	}
	interval subject_range() const
	{
		return interval(j, j + len);
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
	DiagonalSeeds intersect(const DiagonalSeeds &x) const
	{
		if (diag() != x.diag())
			return DiagonalSeeds();
		else {
			const interval q = ::intersect(query_range(), x.query_range());
			return DiagonalSeeds(q.begin_, ::intersect(subject_range(), x.subject_range()).begin_, q.length(), 0);
		}
	}
	bool is_enveloped(const DiagonalSeeds &x) const
	{
		return score <= x.score
			&& query_range().overlap_factor(x.query_range()) == 1
			&& subject_range().overlap_factor(x.subject_range()) == 1;
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
	static bool cmp_subject_end(const DiagonalSeeds &x, const DiagonalSeeds &y)
	{
		return x.subject_end() < y.subject_end();
	}
	static bool cmp_heuristic(const DiagonalSeeds &x, const DiagonalSeeds &y)
	{
		return (x.subject_end() < y.subject_end() && x.j < y.j)
			|| (x.j - y.j < y.subject_end() - x.subject_end());
	}
	static bool cmp_diag(const DiagonalSeeds &x, const DiagonalSeeds &y)
	{
		return x.diag() < y.diag() || (x.diag() == y.diag() && x.j < y.j);
	}
	friend int abs_shift(const DiagonalSeeds &x, const DiagonalSeeds &y)
	{
		return abs(x.diag() - y.diag());
	}
	friend std::ostream& operator<<(std::ostream &s, const DiagonalSeeds &d)
	{
		s << "i=" << d.i << " j=" << d.j << " l=" << d.len << " score=" << d.score;
		return s;
	}
	int i, j, len, score;//query_pos, subject_pos
	hit hit_;
	string qry_;
	string sbj_;
};
/*
struct DiagonalSegment
{
	DiagonalSegment():
		len(0)
	{}

	DiagonalSegment(const TranslatedPosition &i, int j, int len, int score=0):
		i(i),
		j(j),
		len(len),
		score(score)
	{}

	DiagonalSegment(const DiagonalSegment &d, Frame frame):
		i(TranslatedPosition(d.i, frame)),
		j(d.j),
		len(d.len),
		score(d.score)
	{}

	int subject_last() const
	{
		return j + len - 1;
	}

	TranslatedPosition query_last() const
	{
		return i + len - 1;
	}

	int subject_end() const
	{
		return j + len;
	}

	TranslatedPosition query_end() const
	{
		return i + len;
	}

	int diag() const
	{
		return i - j;
	}

	DiagonalSegment& set_score(const TranslatedSequence &query, const sequence &subject)
	{
		TranslatedPosition i1 = i;
		int j1 = j;
		for (; j1 < subject_end(); ++j1)
			score += score_matrix(query[i1], subject[j1]);
		return *this;
	}


	friend std::ostream& operator<<(std::ostream &s, const DiagonalSegment &d)
	{
		s << "i=(" << d.i << ") j=" << d.j << " len=" << d.len << " score=" << d.score << std::endl;
		return s;
	}

	interval query_absolute_range(int dna_len) const
	{
		return TranslatedPosition::absolute_interval(i, i + len, dna_len);
	}

	interval query_in_strand_range() const
	{
		return interval(i.in_strand(), (i + len).in_strand());
	}

	interval subject_range() const
	{
		return interval(j, j + len);
	}

	int partial_score(const DiagonalSegment &d) const
	{
		const double overlap = std::max(subject_range().overlap_factor(d.subject_range()), query_in_strand_range().overlap_factor(d.query_in_strand_range()));
		return int((1.0 - overlap)*score);
	}

	void cut_out(const DiagonalSegment &d)
	{
		const int ll = std::min(d.i.translated - i.translated, d.j - j),
			lr = std::min(query_end().translated - d.query_end().translated, subject_end() - d.subject_end());
		int len2;
		if (ll > 0 && ll >= lr) {
			len2 = std::min(len, ll);
		}
		else if (lr > 0 && lr >= ll) {
			len2 = std::min(len, lr);
			i = query_end() - len2;
			j = subject_end() - len2;
		}
		else {
			len2 = 0;
		}
		score = int((double)len2 / len * score);
		len = len2;
	}

	TranslatedPosition i;
	int j, len, score;
};

*/


#endif // __DIAGONALSEGMENT_H__