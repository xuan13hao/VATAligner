#ifndef __SPLICEDSEED_H__
#define __SPLICEDSEED_H__

#include <vector>
#include <algorithm> 
#include <string>
#include <math.h>
#include <assert.h>
#include <set>
#include "../basic/Seed.h"
using std::vector;
/**
 * 	
 * int i, j, len, score;//query_pos, subject_pos
	string qry_;
	string sbj_;
	int qry_id, sbj_id;


template<typename _locr, typename _locl>
vector<DiagonalSeeds<_locr,_locl> > findSpliceSeeds(vector<DiagonalSeeds<_locr,_locl> > & diagonal_segment,int max_gap)
{
    int MAX = 999999;
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
            int score = diagonal_segment[j].len- gap_cost + splice_score;
            if (pressor_score < score)
            {
                pressor_score = score;
                max_pre = j;
            }
        }
        dp[i+1] = pressor_score;
        pre[i+1] = max_pre;
    }
    for (int i = s_; i > 0; ) 
    {
        if (pre[i] != -1) {
            chained_seed.push_back(diagonal_segment[pre[i]]);
            i = pre[i];
        } else {
            i--;
        }
    }
    std::reverse(chained_seed.begin(),chained_seed.end());
    // cout<<"chained_seed = "<<chained_seed.size()<<endl;
    return chained_seed;
}*/
/*
The findSpliceSeeds function takes a vector of seeds and the maximum gap allowed between seeds to form chains. 
It iterates over the seeds in reverse order, computing the maximum right endpoint for each seed by comparing the gaps between consecutive seeds. 
*/

template<typename _locr, typename _locl>
std::vector<DiagonalSeeds<_locr, _locl>> findSpliceSeeds(std::vector<DiagonalSeeds<_locr, _locl>>& seeds, int maxGap)
{
    std::vector<int> dp(seeds.size());
    std::vector<int> prev(seeds.size(), -1);
    std::vector<int> maxLen(seeds.size());
    std::vector<int> maxIdx(seeds.size());
    std::sort(seeds.begin(),seeds.end(),DiagonalSeeds<_locr,_locl>::cmp_subject_end);
    int MAX = 999999;
    int bestScore = 0;
    int bestIdx = -1;
    for (int i = 0; i < seeds.size(); ++i) {
        dp[i] = seeds[i].len;
        maxLen[i] = seeds[i].len;
        maxIdx[i] = i;

        for (int j = 0; j < i; ++j) {
            int gap = seeds[i].j - seeds[j].j - seeds[j].len;
            if (gap > maxGap)
                continue;
            int splice_score = IsSpliceJunction(seeds[j]) ? -MAX : 0;
            int circ_score = 0;
            if(VATParameters::circ)
            {
                // cout<<"circ"<<endl;
               circ_score =  isCandidateBackSplicedJunction(seeds[j])? -MAX : 0;
            }
            int score = dp[j] + seeds[i].len+splice_score+circ_score;
            if (score > dp[i]) {
                dp[i] = score;
                prev[i] = j;
            }
        }

        if (dp[i] > bestScore) {
            bestScore = dp[i];
            bestIdx = i;
        }
    }

    std::vector<DiagonalSeeds<_locr, _locl>> chainedSeeds;
    while (bestIdx >= 0) {
        chainedSeeds.push_back(seeds[bestIdx]);
        bestIdx = prev[bestIdx];
    }

    std::reverse(chainedSeeds.begin(), chainedSeeds.end());
    return chainedSeeds;
}

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
    int gaps = q_distance - s_distance;
    if(gaps == 0)
    {
        return 0;
    }else if (gaps > 150)
    {
        return 1000000;
    }else
    {   //18 means the mean seed length
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
bool IsSpliceJunction(const DiagonalSeeds<_locr,_locl> &i) 
{
    // int i_len = i.len,j_len = j.len;
    string ref_seed_i = i.sbj_str;
    // string ref_seed_j = j.sbj_str;
    bool i_spjunction = false;
    // bool j_spjunction = false;
    if (ref_seed_i.substr(0, 2) == "GT" && ref_seed_i.substr(ref_seed_i.length() - 2, 2) == "AG") {
            i_spjunction=  true;
    }
    // if (ref_seed_j.substr(0, 2) == "GT" && ref_seed_j.substr(ref_seed_j.length() - 2, 2) == "AG") {
    //         j_spjunction=  true;
    // }
    if(i_spjunction)
    {
        return true;
    }
    else
    {
        return false;
    }
}
template<typename _locr, typename _locl>
bool isCandidateBackSplicedJunction(const DiagonalSeeds<_locr,_locl> &j) {
    // Check seed length
    if (j.len< 8) {
        return false; // Seed is too short to be a candidate junction
    }
    // string ref_seed_i = i.sbj_str;
    string qry_seed_j = j.qry_str;
    // Check for splice site motifs
    size_t donorPos = qry_seed_j.find("GT");
    size_t acceptorPos = qry_seed_j.find("AG");

    // Check if both donor and acceptor motifs are present
    if (donorPos == std::string::npos || acceptorPos == std::string::npos) {
        return false; // Seed does not contain both donor and acceptor motifs
    }

    // Check orientation
    if (donorPos > acceptorPos) {
        return true; // Seed has the reverse orientation of splice sites
    }

    return false; // Seed does not meet the criteria for a candidate junction
}

vector<pair<int, int>> find_exons(string sequence) {
    vector<pair<int, int>> exon_regions; // stores start and end positions of exon regions
    int seq_length = sequence.length();
    bool in_exon = false;
    int exon_start = 0;

    for (int i = 0; i < seq_length; i++) {
        if (sequence[i] == 'G' && sequence[i+1] == 'T' && sequence[i+2] == 'A' && sequence[i+3] == 'G') {
            // stop codon, end of exon
            if (in_exon) {
                exon_regions.push_back(make_pair(exon_start, i+3));
                in_exon = false;
            }
        } else if (sequence[i] == 'A' && sequence[i+1] == 'T' && sequence[i+2] == 'G') {
            // start codon, beginning of exon
            in_exon = true;
            exon_start = i;
        }
    }

    if (in_exon) {
        // sequence ends in exon
        exon_regions.push_back(make_pair(exon_start, seq_length-1));
    }

    return exon_regions;
}

#endif // __SPLICEDSEED_H__