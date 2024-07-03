#include <iostream>
#include <vector>
#include <algorithm>
// #include "../../basic/Seed.h"
using std::vector;
// template<typename _locr, typename _locl>
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
		// hit_ = rhs.hit_;
		// qry_id = rhs.qry_id;
		// sbj_id = rhs.sbj_id;
		score = rhs.score;
		len = rhs.len;
		// qry_str = rhs.qry_str;
		// sbj_str = rhs.sbj_str;

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
	friend std::ostream& operator<<(std::ostream &s, const DiagonalSeeds &d)
	{
		s << "i=" << d.i << " j=" << d.j << " l=" << d.len << " score=" << d.score;
		return s;
	}
	int i, j, len, score;//query_pos, subject_pos
	// string qry_str;
	// string sbj_str;
	// int qry_id, sbj_id;
};


vector<DiagonalSeeds> findWholeGenSeeds(vector<DiagonalSeeds>& diagonal_segment, int max_gap)
{
    vector<DiagonalSeeds> chained_seed;
    const int s_ = diagonal_segment.size();
    vector<int> dp(s_+1, 0);
    vector<int> pre(s_+1, -1);
    dp[0] = diagonal_segment[0].score;
    pre[0] = -1;
    sort(diagonal_segment.begin(), diagonal_segment.end(), DiagonalSeeds::cmp_subject_end);

    for (auto i = 1; i < s_; i++)
    {
        int pressor_score = diagonal_segment[i].score;
        int max_pre = i;
        const int subject_end = diagonal_segment[i].j;
        for (auto j = i - 1; j >= 0 && i - j <= 50; j--)
        {
            const int diff = subject_end - diagonal_segment[j].j;
            if (diff > max_gap) break;
            const int score_j = diagonal_segment[j].score;
            if (pressor_score < score_j)
            {
                pressor_score = score_j;
                max_pre = j;
            }
        }
        dp[i] = pressor_score;
        pre[i] = max_pre;
    }
    for (int i = s_; i > 0; ) {
        if (pre[i] != -1) {
            chained_seed.push_back(diagonal_segment[pre[i]]);
            i = pre[i];
        } else {
            i--;
        }
    }
    reverse(chained_seed.begin(), chained_seed.end());
    return chained_seed;
}

int main() {
    std::vector<DiagonalSeeds> diagonal_segment = {
        { 0, 0, 5, 3 },
        { 1, 2, 6, 5 },
        { 3, 6, 7, 4 },
        { 4, 9, 4, 2 },
        { 7, 13, 5, 4 }
    };
    int max_gap = 3;

    std::vector<DiagonalSeeds> chained_seeds = findWholeGenSeeds(diagonal_segment, max_gap);

    for (const auto& seed : chained_seeds) {
        std::cout << "i: " << seed.i << ", j: " << seed.j << ", len: " << seed.len << ", score: " << seed.score << std::endl;
    }

    return 0;
}