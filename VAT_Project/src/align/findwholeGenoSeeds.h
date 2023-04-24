#ifndef __FINDWHOLEGENOSEEDS_H__
#define __FINDWHOLEGENOSEEDS_H__

#include <vector>
#include <algorithm> 
#include <string>
#include <cmath>
#include <assert.h>
#include <set>
#include "../basic/Seed.h"
#define MAX 999999999
// #include "kvec.h"
using std::vector;
/**
 * 	
 * int i, j, len, score;//query_pos, subject_pos
	string qry_;
	string sbj_;
	int qry_id, sbj_id;


vector<int> findPossibleExons(string dna, vector<int>& seeds, int k) {
    int n = dna.size();
    vector<int> dp(n+1, 0);
    vector<int> pre(n+1, -1);

    for (int i = 0; i < n; i++) {
        if (i < k-1) {
            dp[i+1] = 0; // non-coding region
        } else {
            dp[i+1] = (dna[i] == 'A' || dna[i] == 'T' || dna[i] == 'C' || dna[i] == 'G') ? 1 : 0; // coding region
        }
    }

    for (int i = 1; i <= n; i++) {
        for (int j = 0; j < i; j++) {
            if (find(seeds.begin(), seeds.end(), j+1) != seeds.end()) { // check if seed is qualified
                if (dp[j+1] > 0 && dp[j+1] + 1 > dp[i+1]) {
                    dp[i+1] = dp[j+1] + 1; // update max exon length
                    pre[i+1] = j; // mark previous position
                }
            }
        }
    }

    vector<int> exons;
    for (int i = n; i > 0; ) {
        if (pre[i] != -1) {
            exons.push_back(pre[i]);
            i = pre[i];
        } else {
            i--;
        }
    }
    reverse(exons.begin(), exons.end());
    return exons;
}
*/

template<typename _locr, typename _locl>
vector<DiagonalSeeds<_locr,_locl> > findWholeGenSeeds(vector<DiagonalSeeds<_locr,_locl> > & diagonal_segment,int max_gap)
{
    vector<DiagonalSeeds<_locr,_locl> > chained_seed;
    int s_ =  diagonal_segment.size();
    int i, score, qdiff, tdiff, diffdiff, gap_cost, j, best_j;
    int h = 50; // number of previous anchors to check
    int match_score = 0;
    vector<int> dp(s_, 0);
    vector<int> pre(s_, -1);
    dp[0] = diagonal_segment[0].score;
    std::sort(diagonal_segment.begin(),diagonal_segment.end(),DiagonalSeeds<_locr,_locl>::cmp_subject_end);
    
    for (size_t i = 1; i < diagonal_segment.size(); i++)
    {
        int pressor_score = diagonal_segment[i].score; 
        int max_pre = i;
        // int gap_score = -MAX;
        for(j = 0; j < i; j++) 
        {
            // int gap_score = -MAX;
            // int gap_distance = distanceSegment(diagonal_segment[j], diagonal_segment[i]);
                // if it fits the requirement of macro-gap
            // int tdiff = static_cast<int>(diagonal_segment[j].j) - static_cast<int>(diagonal_segment[i].j);
            // int qdiff = static_cast<int>(diagonal_segment[j].i) - static_cast<int>(diagonal_segment[i].i) ;
            // if(abs(tdiff - qdiff) <= max_gap)
            // {
            //     gap_score = diagonal_segment[i].score_; 
            // }
            // int gap_cost = gapCost(diagonal_segment[j], diagonal_segment[i]);    
            // int splice_score = IsSpliceJunction(diagonal_segment[j], diagonal_segment[i]) ? -MAX : 0; 
            // int score = diagonal_segment[j].len+ gap_distance - gap_cost + splice_score;
            int score = diagonal_segment[j].score - gap_cost
            if (pressor_score < score)
            {  
                pressor_score = score;
                max_pre = j;
            } 
        }
        dp[i] = pressor_score;
        pre[i] = max_pre;
    }
    vector<int> exons;
    for (int i = s_; i > 0; ) {
        if (pre[i] != -1) {
            chained_seed.push_back(diagonal_segment[pre[i]]);
            exons.push_back(pre[i]);
            i = pre[i];
        } else {
            i--;
        }
    }
    std::reverse(exons.begin(),exons.end());
    std::reverse(chained_seed.begin(),chained_seed.end());

    return chained_seed;

}
/**
 * int i, j, len, score;//query_pos, subject_pos
*/
template<typename _locr, typename _locl>
int distanceSegment(const DiagonalSeeds<_locr,_locl> &j, const DiagonalSeeds<_locr,_locl> &i) //alpha(j,i)
{
    int q_distance = i.i - j.i;
    int s_distance = i.j - j.j;
    int distance = min(q_distance,s_distance);
    int i_len = i.len;
    int alph_value = min(distance,i_len);
    return alph_value;
}

template<typename _locr, typename _locl>
int gapCost(const DiagonalSeeds<_locr,_locl> &j, const DiagonalSeeds<_locr,_locl> &i) //beta(j,i)
{
    int q_distance = i.i - j.i;
    int s_distance = i.j - j.j;
    int gaps = abs(q_distance - s_distance);
    if(gaps == 0)
    {
        return 0;
    }else if (gaps > 150)
    {
        return 1000000;
    }else
    {   //15 means the mean seed length
        int gapcost = 0.01*15*abs(gaps)+0.05*log2(abs(gaps));
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

template<typename _locr, typename _locl>
bool IsSpliceJunction(const DiagonalSeeds<_locr,_locl> &j, const DiagonalSeeds<_locr,_locl> &i) 
{
    string ref_seed_i, ref_seed_j, query_i,query_j;
    ref_seed_i = i.sbj_str;
    ref_seed_j = j.sbj_str;

    query_i = i.qry_str;
    query_j = j.qry_str;

    bool ref_i_spjunction = false;
    bool ref_j_spjunction = false;
    bool query_i_spjunction = false;
    bool query_j_spjunction = false;

    if (ref_seed_i.substr(0, 2) == "GT" && ref_seed_i.substr(ref_seed_i.length() - 2, 2) == "AG") {
            bool ref_i_spjunction=  true;
    }
    if (ref_seed_j.substr(0, 2) == "GT" && ref_seed_j.substr(ref_seed_j.length() - 2, 2) == "AG") {
            bool ref_j_spjunction=  true;
    }
    if (query_i.substr(0, 2) == "GT" && query_i.substr(query_i.length() - 2, 2) == "AG") {
            bool query_i_spjunction=  true;
    }
    if (query_j.substr(0, 2) == "GT" && query_j.substr(query_j.length() - 2, 2) == "AG") {
            bool query_j_spjunction=  true;
    }

    if((query_i_spjunction && ref_i_spjunction)  || (ref_j_spjunction&&query_j_spjunction))
    {
        return true;
    }
    else
    {
        return false;
    }

}


#endif // __FINDWHOLEGENOSEEDS_H__