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
/**
 * 	
 * int i, j, len, score;//query_pos, subject_pos
	string qry_;
	string sbj_;
	int qry_id, sbj_id;
*/

template<typename _locr, typename _locl>
vector<DiagonalSeeds<_locr,_locl> > findSpliceSeeds(vector<DiagonalSeeds<_locr,_locl> > & diagonal_segment,int max_gap)
{
    vector<DiagonalSeeds<_locr,_locl> > chained_seed;
    int s_ =  diagonal_segment.size();
    int i, score, qdiff, tdiff, diffdiff, gap_cost, j, best_j;
    int h = 50; // number of previous anchors to check
    int match_score = 0;
    int *p1_score = new int [s_]; 
    int *p1_track = new int [s_];

    p1_score[0] = diagonal_segment[0].len
    p1_track[0] = -1;
    
    std::sort(diagonal_segment.begin(),diagonal_segment.end(),DiagonalSeeds<_locr,_locl>::cmp_subject_end);
    for (size_t i = 1; i < diagonal_segment.size(); i++)
    {
        // cout<<"findSpliceSegments"<<endl;
        int pressor_score = diagonal_segment[i].len; 
        int max_pre = i;
        // int gap_score = -MAX;
        for(j = 0; j < i; j++) 
        {
            // int gap_score = -MAX;
            int gap_distance = distanceSegment(diagonal_segment[j], diagonal_segment[i]);
                // if it fits the requirement of macro-gap
            //int ref_gap = static_cast<int>(diagonal_segment[j].j) - static_cast<int>(diagonal_segment[i].j) - static_cast<int>(diagonal_segment[i].len);
            //int qry_gap = static_cast<int>(diagonal_segment[j].i) - static_cast<int>(diagonal_segment[i].i) - static_cast<int>(diagonal_segment[i].len);
            int tdiff = static_cast<int>(diagonal_segment[j].j) - static_cast<int>(diagonal_segment[i].j);
            int qdiff = static_cast<int>(diagonal_segment[j].i) - static_cast<int>(diagonal_segment[i].i) ;
            // if(abs(tdiff - qdiff) <= max_gap)
            // {
            //     gap_score = diagonal_segment[i].score_; 
            // }
            int gap_cost = gapCost(diagonal_segment[j], diagonal_segment[i]);    
            int splice_score = IsSpliceJunction(diagonal_segment[j], diagonal_segment[i]) ? -MAX : 0; 
            int score = diagonal_segment[j].len+ gap_distance - gap_cost + splice_score;
            if (pressor_score < score)
            {  
                pressor_score = score;
                max_pre = j;
            } 
        }
        p1_score[i] = pressor_score;
        p1_track[i] = max_pre;
    }
    set<int> visited;
    set<int>::iterator it;
    for (size_t i = 0; i < s_; i++)
    {
        visited.insert(p1_track[i]);
    }
    for (it = visited.begin(); it != visited.end(); ++it)
    {
        // cout<<*it<<endl;
        if (*it == -1)
        {
            chained_seed.push_back(diagonal_segment[0]);
        }else
        {
            chained_seed.push_back(diagonal_segment[*it]);
        }   
        
    }

    // std::reverse(chained_seed.begin(),chained_seed.end());
    delete [] p1_score;
    delete [] p1_track;

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
bool IsSpliceJunction(const DiagonalSeeds<_locr,_locl> &j, const DiagonalSeeds<_locr,_locl> &i) 
{
    const SequenceSet<_val> *ref = ReferenceSeqs<_val>::data_;
    const _val* sbj_j = ref->data(j.i);
    const _val* sbj_i = ref->data(i.i);
    int i_n = 0,j_n = 0;;
    int i_len = i.len,j_len = j.len;
    string ref_seed_i, ref_seed_j;
    string splice_signal_AG("AG");
    string splice_signal_GT("GT");
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
    if ((ref_seed_i.find(splice_signal_AG) != std::string::npos 
    && ref_seed_j.find(splice_signal_GT) != std::string::npos)
    )
    {
        return true;
    }
    // delete ref1;
    return false;
}


#endif // __SPLICEDSEED_H__