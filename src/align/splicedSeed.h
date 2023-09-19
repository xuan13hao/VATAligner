#ifndef __SPLICEDSEED_H__
#define __SPLICEDSEED_H__

#include <vector>
#include <algorithm> 
#include <string>
#include <math.h>
#include <assert.h>
#include <set>
#include "../basic/Seed.h"
#define MAX 999999999
// #include "kvec.h"
using std::vector;


template<typename _val>
vector<Segment<_val> > findSpliceSegments(vector<Segment<_val> > & diagonal_segment,int max_gap)
{
    vector<DiagonalSeeds<_locr,_locl> > chained_seed;
    int s_ =  diagonal_segment.size();
    vector<int> dp(s_+1, 0);
    vector<int> pre(s_+1, -1);
    dp[0] = diagonal_segment[0].len;
    std::sort(diagonal_segment.begin(),diagonal_segment.end(),DiagonalSeeds<_locr,_locl>::cmp_subject_end);
    
    for (auto i = 1; i < s_; i++)
    {
        int pressor_score = diagonal_segment[i].len;
        int max_pre = i;
        const int subject_end = diagonal_segment[i].j;
        for (auto j = i - 1; j >= 0 && i - j <= 50; j--)
        {
            const int diff = subject_end - diagonal_segment[j].j;
            if (diff > max_gap) break;
            
            int gap_cost = gapCost(diagonal_segment[j], diagonal_segment[i]);    
            int splice_score = IsSpliceJunction(diagonal_segment[j], diagonal_segment[i]) ? -MAX : 0; 
            int score = diagonal_segment[j].traceback_->len_+ gap_distance - gap_cost + splice_score;
            if (pressor_score < score)
            {
                pressor_score = score;
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
    std::reverse(chained_seed.begin(),chained_seed.end());
    return chained_seed;
}
template<typename _val>
int distanceSegment(const Segment<_val> &j, const Segment<_val> &i) //alpha(j,i)
{
    int q_distance = i.traceback_->query_begin_ - j.traceback_->query_begin_;
    int s_distance = i.traceback_->subject_begin_ - j.traceback_->subject_begin_;
    int distance = min(q_distance,s_distance);
    int i_len = i.traceback_->len_;
    int alph_value = min(distance,i_len);
    return alph_value;
}

template<typename _val>
int gapCost(const Segment<_val> &j, const Segment<_val> &i) //beta(j,i)
{
    int q_distance = i.traceback_->query_begin_ - j.traceback_->query_begin_;
    int s_distance = i.traceback_->subject_begin_ - j.traceback_->subject_begin_; 
    int gaps = q_distance - s_distance;
    if(gaps == 0)
    {
        return 0;
    }else if (gaps > 150)
    {
        return 1000000;
    }else
    {   //18 means the mean seed length
        int gapcost = 0.01*18*abs(gaps)+0.05*log2(abs(gaps));
        return gapcost;
    }
    
}
/**
Use std::string::find as follows:

if (s1.find(s2) != std::string::npos) {
    std::cout << "found!" << '\n';
}
Note: "found!" will be printed if s2 is a substring of s1, both s1 and s2 are of type std::string.

splice juntion GT.........AG 
*/
template<typename _val>
bool IsSpliceJunction(const Segment<_val> &j, const Segment<_val> &i) 
{
    const SequenceSet<_val> *ref = ReferenceSeqs<_val>::data_;
    const _val* sbj_j = ref->data(j.traceback_->subject_begin_);
    const _val* sbj_i = ref->data(i.traceback_->subject_begin_);
    int i_n = 0,j_n = 0;;
    int i_len = i.traceback_->len_,j_len = j.traceback_->len_;
    string ref_seed_i, ref_seed_j;
    // string splice_signal_AG("AG");
    // string splice_signal_GT("GT");
    while(*sbj_i != AlphabetSet<_val>::PADDING_CHAR && i_n < i_len) 
    {
		++sbj_i;
        ref_seed_i.push_back(AlphabetAttributes<_val>::ALPHABET[*(++sbj_i)]);
		++i_n;
	}
    while(*sbj_j != AlphabetSet<_val>::PADDING_CHAR && j_n < j_len) 
    {
		++sbj_j;
        ref_seed_j.push_back(AlphabetAttributes<_val>::ALPHABET[*(++sbj_j)]);
		++j_n;
	}
    bool i_spjunction = false;
    bool j_spjunction = false;
    if (ref_seed_i.substr(0, 2) == "GT" && ref_seed_i.substr(ref_seed_i.length() - 2, 2) == "AG") {
            bool i_spjunction=  true;
    }
    if (ref_seed_j.substr(0, 2) == "GT" && ref_seed_j.substr(ref_seed_j.length() - 2, 2) == "AG") {
            bool j_spjunction=  true;
    }
    if(i_spjunction || j_spjunction)
    {
        return true;
    }
    else
    {
        return false;
    }
}
#endif // __SPLICEDSEED_H__