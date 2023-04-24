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
	// DiagonalSeeds(int query_pos, int subject_pos, int len, int score, Hits<_locr,_locl>& h, int& q_, int& s_,string& q, string& s) :
	// i(query_pos),
	// j(subject_pos),
	// len(len),
	// score(score),
	// hit_(h),
	// qry_id(s_),
	// sbj_id(q_),
	// qry_str(q),
	// sbj_str(s)
	// {
		
	// }

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


// vector<DiagonalSeeds> findWholeGenSeeds(vector<DiagonalSeeds>& diagonal_segment, int max_gap)
// {
//     vector<DiagonalSeeds> seeds;
//     int n = diagonal_segment.size();
//     int i = 0;

//     while (i < n) {
//         // Start a new seed
//         int j = i + 1;
//         int start = diagonal_segment[i].i;
//         int end = diagonal_segment[i].i + diagonal_segment[i].len - 1;

//         while (j < n && diagonal_segment[j].i <= end + max_gap) {
//             // Extend the seed
//             end = diagonal_segment[j].i + diagonal_segment[j].len - 1;
//             j++;
//         }

//         // Add the seed to the list
//         seeds.push_back(DiagonalSeeds(start, diagonal_segment[i].j, end - start + 1, 0));

//         // Move to the next diagonal
//         i = j;
//     }

//     return seeds;
// }

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
        dp[i+1] = pressor_score;
        pre[i+1] = max_pre;
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

int main()
{
    // Create some diagonal segments
    vector<DiagonalSeeds> diagonal_segments = {
        DiagonalSeeds(10, 20, 30, 50),
        DiagonalSeeds(40, 50, 10, 30),
        DiagonalSeeds(70, 80, 20, 10),
        DiagonalSeeds(110, 120, 15, 40),
        DiagonalSeeds(140, 150, 25, 20),
        DiagonalSeeds(190, 200, 30, 80)
    };

    // Find whole genome seeds
    vector<DiagonalSeeds> whole_genome_seeds = findWholeGenSeeds(diagonal_segments, 5);

    // Print the seeds
    for (const auto& seed : whole_genome_seeds) {
        std::cout << "Seed: i=" << seed.i << " j=" << seed.j << " len=" << seed.len << " score=" << seed.score << std::endl;
    }

    return 0;
}
