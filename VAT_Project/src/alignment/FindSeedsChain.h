#ifndef __FINDSEEDSCHAIN_H__
#define __FINDSEEDSCHAIN_H__
#define MAX 65535
#include <vector>
#include <algorithm> 
#include <string.h>
#include <math.h>
#include <assert.h>
#include "../basic/DiagonalSeeds.h"

// #include "kvec.h"
using std::vector;
// defines a chain of colinear seeds
// chain contains the set of seeds, is_splice indicates whether the two adjacent seeds correspond to a splicing site
struct SeedChainType   {
    std::vector<DiagonalSeeds> chain;
    std::vector<bool> is_splice;
    int num_seeds;
    int score;

    SeedChainType& operator =(const SeedChainType& a)
    {
        chain = a.chain;
        is_splice = a.is_splice;
        num_seeds = a.num_seeds;
        score = a.score;
        return *this;
    }
};





vector<DiagonalSeeds> findOptimalSeeds(vector<DiagonalSeeds>& diagonal_segment,int q_len,int max_gap)
{
    vector<DiagonalSeeds> chained_seed;
    // uint32_t target = 0;
    int s_ =  diagonal_segment.size();
    // if(s_ == 0)
    // {
    //     return chained_seed;
    // }
    // if (s_ == 1)
    // {
    //     chained_seed.push_back(diagonal_segment[0]);
    // }
    int i, score, qdiff, tdiff, diffdiff, gap_cost, j, best_j;
    int ref_offset = 0;
    int h = 50; // number of previous anchors to check
    // int ref_offset = 0; // offset into the scores vector of the current target - sets the backward limit for finding chained anchors
    int match_score = 0;
    int pre[s_];
    int f[s_];
    int track[s_];
    f[0] = diagonal_segment[0].len;
    int *p1_score = new int [s_]; 
    int *p1_track = new int [s_];
    p1_score[0] = diagonal_segment[0].score;
    p1_track[0] = -1;
    // std::sort(diagonal_segment.begin(),diagonal_segment.end(),DiagonalSeeds::cmp_subject_end);
    for (size_t i = 1; i < diagonal_segment.size(); i++)
    {
        // int weight = diagonal_segment[i].len;
        int tmp = diagonal_segment[i].score; 
        track[0] = 0;
        int max_pre = i;
        for(j = 0; j < i; j++) 
        {
            int tdiff = static_cast<int>(diagonal_segment[j].j) - static_cast<int>(diagonal_segment[i].j) ;
            int qdiff = static_cast<int>(diagonal_segment[j].i) - static_cast<int>(diagonal_segment[i].i) ;
            int ref_gap = static_cast<int>(diagonal_segment[j].j) - static_cast<int>(diagonal_segment[i].j) - static_cast<int>(diagonal_segment[i].len);
            int qry_gap = static_cast<int>(diagonal_segment[j].i) - static_cast<int>(diagonal_segment[i].i) - static_cast<int>(diagonal_segment[i].len);

            if(tdiff <= 0 || qdiff <= 0 || qdiff > max_gap || tdiff > max_gap) 
            { // we often have multiple hits to the same target pos, so we can't let them chain with either other
                continue;
            }
            // if(ref_gap <= 0 || qry_gap <= 0 || qry_gap > max_gap || ref_gap > max_gap) 
            // { // we often have multiple hits to the same target pos, so we can't let them chain with either other
            //         continue;
            // }
            match_score = diagonal_segment[j].score;
            gap_cost = (diffdiff == 0 ? 0 : 0.01 * match_score * diffdiff + 0.5 * log2(diffdiff)); // affine gap penalty a la minimap2
            qdiff = (qdiff > tdiff ? tdiff : qdiff); // now not qdiff but the minimum difference
            score = diagonal_segment[j].score + (qdiff > match_score ? match_score : qdiff) - gap_cost; 
            
        //    score = diagonal_segment[j].score;
             if (tmp <= score)
             {  
                 tmp = score;
                 max_pre = j;
             } 
        }
        p1_score[i] = tmp;
        p1_track[i] = max_pre;
    }
    
    for (size_t i = 0; i < s_; i++)
    {
        cout<<p1_track[i]<<"\t"<<p1_score[i]<<endl;
    }
    cout<<endl;
    bool *tracked = new bool [s_];
    memset(tracked, false, s_ * sizeof(bool));
    cout<<"s = "<<s_<<endl;
    
    bool flag = true;
    for (size_t j = s_ -1; j >= 0; j--)
    { 
        
        if (tracked[j])
        {
            continue;
        }
        int c = j;
        if (p1_track[c] == -1)
        {
            chained_seed.push_back(diagonal_segment[c]);
            break;
        }
        while (true)
        {
            tracked[c] = true;
            chained_seed.push_back(diagonal_segment[c]);
            if (p1_track[c] == c)
            {
                break;
            }
            else
            {
                c = p1_track[c];
            } 
        }
    }
    std::reverse(chained_seed.begin(),chained_seed.end());
    cout<<"chain seed = "<<chained_seed.size()<<endl;
    // cout<<"here"<<endl;
    delete [] tracked;
    delete [] p1_score;
    delete [] p1_track;

    // int k = max_element(f,f+diagonal_segment.size())-f;
    
    // bool flag = true;
    // int len = f[k];
    // cout<<k<<"\t"<<len<<endl;
    // while (flag)
    // {
    //     chained_seed.push_back(diagonal_segment[k]);
    //     k = pre[k];
    //     if (k == 0)
    //     {
    //         flag = false;
    //     }
    //     len --;
    // }
    // cout<<"chain seed = "<<chained_seed.size()<<endl;
//    delete [] tracked;
    return chained_seed;
}


void findSeedChain(vector<DiagonalSeeds>& diagonal_segment,vector<SeedChainType>& seed_chains,int q_len)
{
    // vector<SeedChainType> seed_chains;
    if (diagonal_segment.size() <=0)
    {
        return;
    }
    int qry_map_len = 0.9*q_len;
    for (size_t i = 0; i < diagonal_segment.size(); i++)
    {
        // if (diagonal_segment[i].len >= qry_map_len)
        // {
            SeedChainType sc;
            sc.chain.push_back(diagonal_segment[i]);
            sc.num_seeds = 1;
            seed_chains.push_back(sc);
        // }
    }
    if (seed_chains.size() > 0)
    {
       return;
    }
    std::sort(diagonal_segment.begin(),diagonal_segment.end(),DiagonalSeeds::cmp_subject_end);
    // construct the pre-set for chaining, which contains the set of pre-seeds for each seed
    int **pre_seed_id = new int *[diagonal_segment.size()];   
    int **pre_gap_score = new int *[diagonal_segment.size()];      // assuming the seed break is a gap
    int **pre_spl_score = new int *[diagonal_segment.size()];      // assuming the seed break is a splice site
    for(int i = 0; i < diagonal_segment.size(); ++ i) {
        pre_seed_id[i] = new int [diagonal_segment.size()];
        pre_gap_score[i] = new int [diagonal_segment.size()];
        pre_spl_score[i] = new int [diagonal_segment.size()];
    }
    int *pre_count = new int [diagonal_segment.size()];    // stores the number of pre-seeds for each seed
    memset(pre_count, 0, diagonal_segment.size() * sizeof(int));
    int seed_break, gap_size, overlap;
    for(int i = 0; i < diagonal_segment.size() - 1; ++ i)   {
        int bound = q_len;
        int j = i + 1;
        while(j < diagonal_segment.size()) 
        {
            // set the right bound as the end position of the next seed
            if(diagonal_segment[i].i + diagonal_segment[i].len - diagonal_segment[j].i <= overlap)    {
                bound = diagonal_segment[j].i + diagonal_segment[j].len - 1;
                break;
            }
            ++ j;
        }
        // DEBUG
        //cout << "check index:\t" << i << "\t" << j << "\t" << bound << "\t" << seeds[j].qry_pos << "\t" << pre_count[j] << endl;
        
        // the remaining seeds must overlap with the seed for less than overlap, 
        // because all seeds are sorted by their starting locations
        /*
        int macro_gap_size_;    // the maximum allowed INDEL size (i.e. the difference in alignment length)
        int macro_seed_break_;  // the maximum seed break (i.e. the distance between two seeds)
        int macro_overlap_;     // the maximum seed overlap (i.e. the overlap between two seeds if they can be chained)*/
        for(; j < diagonal_segment.size(); ++ j)   
        {
            
            if( diagonal_segment[j].i <= bound)    
            {
                // here we need to consider splicing event and long-range gap
                int gap_score = -MAX;
                // if it fits the requirement of macro-gap
                int ref_gap = static_cast<int>(diagonal_segment[j].j) - static_cast<int>(diagonal_segment[i].j) - static_cast<int>(diagonal_segment[i].len);
                int qry_gap = static_cast<int>(diagonal_segment[j].i) - static_cast<int>(diagonal_segment[i].i) - static_cast<int>(diagonal_segment[i].len);
                // DEBUG
                //cout << "REF Seed POS:\t" << seeds[j].ref_pos << "\t" << seeds[i].ref_pos << "\t" << seeds[i].ref_len << endl;
                //cout << "QRY Seed POS:\t" << seeds[j].qry_pos << "\t" << seeds[i].qry_pos << "\t" << seeds[i].qry_len << endl;
                //cout << "GAP sizes:\t" << ref_gap << "\t" << qry_gap << endl;

                if(qry_gap <= seed_break && abs(ref_gap - qry_gap) <= gap_size)    
                {
                    int band = 2 * (abs(ref_gap - qry_gap) + 1);
                    gap_score = diagonal_segment[i].score;
                    // gap_score = align_obj.AlignGlobalPairwise(
                    //     ref_seq.sequence_[seeds[i].ref_id] + seeds[i].ref_pos + 1, qry_seq.sequence_[seeds[i].qry_id] + seeds[i].qry_pos + 1,
                    //     ref_gap, qry_gap, band, score
                    // );    
                }
                   
                // int splice_score = IsSpliceJunction(ref_seq, qry_seq, seeds[i], seeds[j]) ? 0 : -MAX;
                int splice_score = 0;
                // record the score
                if(gap_score > -MAX || splice_score > -MAX) {
                    // it means that the seed break might make sence, record
                    pre_seed_id[j][pre_count[j]] = i;
                    pre_gap_score[j][pre_count[j]] = gap_score;
                    pre_spl_score[j][pre_count[j]] = splice_score;
                    ++ pre_count[j];
                }

            }
        }
    }

 // use dynamic programming to find the best macro-gap or splice alignment
    int *p1_score = new int [diagonal_segment.size()]; int *p1_track = new int [diagonal_segment.size()]; bool *p1_splice = new bool [diagonal_segment.size()];
    // initialization
    p1_score[0] = diagonal_segment[0].score; p1_track[0] = -1;
    for(int i = 1; i < diagonal_segment.size(); ++ i)   {
        // perform dynamic programming, we need to record whether the seed break is a splicing site
        int max_score = diagonal_segment[i].score; int max_pre = i; bool is_splice = false;
        for(int j = 0; j < pre_count[i]; ++ j)   {
            int gs = pre_gap_score[i][j] + p1_score[pre_seed_id[i][j]];     // prev gap score
            int ss = pre_gap_score[i][j] + p1_score[pre_seed_id[i][j]];     // prev splice score
            if(gs > max_score)    {
                max_score = gs; max_pre = pre_seed_id[i][j]; is_splice = false;
            }
            if(ss > max_score)    {
                max_score = ss; max_pre = pre_seed_id[i][j]; is_splice = true;
            }
        }
        p1_score[i] = max_score;
        p1_track[i] = max_pre;
        p1_splice[i] = is_splice;
    }

    // traceback and recover the seed chains
    bool *tracked = new bool [diagonal_segment.size()];
    memset(tracked, false, diagonal_segment.size() * sizeof(bool));
    for(int j = diagonal_segment.size() - 1; j >= 0; j --)   {
        if(tracked[j]) continue;
        SeedChainType sc;
        int c = j;
        while(true) {
            tracked[c] = true;
            sc.chain.push_back(diagonal_segment[c]);
            sc.is_splice.push_back(p1_splice[c]);
            ++ sc.num_seeds;
            if(p1_track[c] == c)    {
                break;
            }   else
            {
                c = p1_track[c];
            }
        }
        // reverse the orders of the information in the seed chain
        std::reverse(sc.chain.begin(), sc.chain.end());
        std::reverse(sc.is_splice.begin(), sc.is_splice.end());
        seed_chains.push_back(sc);
    }

    // DEBUG
    //Seeding seed_obj;
    //for(int x = 0; x < seed_chains.size(); ++ x)   {
    //    cout << ">>>>>>>>>>>>>>>> Print sorted seeds: >>>>>>>>>>>>>>>>" << endl;
    //    for(int y = 0; y < seed_chains[x].chain.size(); ++ y)   {
    //        seed_obj.PrintSingleSeed(ref_seq, qry_seq, seed_chains[x].chain[y]);
    //    }
    //}

    // collect memory
    delete [] tracked;
    delete [] p1_score;
    delete [] p1_track;
    delete [] p1_splice;
    for(int i = 0; i < diagonal_segment.size(); ++ i)   {
        delete [] pre_seed_id[i];
        delete [] pre_gap_score[i];
        delete [] pre_spl_score[i];
    }
    delete [] pre_seed_id;        // if (diagonal_segment[i].len >= qry_map_len)
        // {
    
    // DEBUG
    //cout << "DONE collectin memory" << endl;
        // if (diagonal_segment[i].len >= qry_map_len)
        // {

    return;
}
/**
 * @brief 

    def get_chains(self, minimium_chaining_score=20, minimum_minimizers_on_chain=2):
        chaining_scores = [self.w]  # Maximal score up to the ith anchor at ith position
        best_predecessors = {0: -1}  # Best predecessor anchors of anchor at key


        # Formula 1 and 2 from Heng Li's Minimap2 paper
        for i in range(1, len(self.anchors)):
            #print("Finding score for i=%d" % i)
            scores = {}  # Scores (value) between chain j (key) and i
            for j in range(i-1, -1, -1):  # All j's lower than i

                score = max(chaining_scores[j] + self.distance_between_anchors(i, j) - self.anchor_gap_cost(i, j), self.w)
                #score = max(chaining_scores[j] + self.anchors_score(i, j), self.w)
                #print("   j=%d. Dist: %.3f. Score: %.2f. Gap cost: %d" % (j, self.distance_between_anchors(i, j), score, self.anchor_gap_cost(i, j)))
                scores[j] = score
                if j < i - 1:
                    # Using heuristic that chaining this far down probably gives a lower score
                    break
            #print("    scores: %s" % str(scores))
            # Find best predecessor as the one with max score
            best_predecessor, best_predecessor_score = max(scores.items(), key=lambda s: s[1] + 0.000001 * s[0])  # Important to prioritize highest key if equal score
            best_predecessor_score = scores[best_predecessor]
            if best_predecessor_score == self.w:
                best_predecessor = -1
            #print("    Best predecessor of %d: %d" % (i, best_predecessor))
            best_predecessors[i] = best_predecessor
            chaining_scores.append(best_predecessor_score)

        # Backtracking
        chains = []
        #print(" == Backtracking ==")
        used_anchors = set()
        for i in range(len(self.anchors)-1, -1, -1):
            #print("Checking anchor %d" % i)
            if i in used_anchors:
                continue
            current_chain = Chain()
            current_anchor = i
            while True:
                used_anchors.add(current_anchor)
                current_chain.anchors.append(self.anchors[current_anchor])
                # Find best predecessor
                best = best_predecessors[current_anchor]
                #print("  Found best %d" % best)
                if best == -1 or best in used_anchors:
                    break
                current_anchor = best
            current_chain.chaining_score = chaining_scores[i]
            if current_chain.chaining_score >= minimium_chaining_score and \
                    len(current_chain.anchors) >= minimum_minimizers_on_chain:
                # Hack for now, we need to reverse order to be compatible with old setup
                current_chain.anchors = current_chain.anchors[::-1]
                chains.append(current_chain)

 
  
 */

vector<DiagonalSeeds> chainingSeeds(vector<DiagonalSeeds>& diagonal_segment,int q_len,int max_gap)
{
    vector<DiagonalSeeds> chained_seed;
    int s_ =  diagonal_segment.size();
    int i, score, qdiff, tdiff, diffdiff, gap_cost, j, best_j;
    int h = 50; // number of previous anchors to check
    int match_score = 0;
    int *p1_score = new int [s_]; 
    int *p1_track = new int [s_];

    p1_score[0] = diagonal_segment[0].score;
    p1_track[0] = -1;

    // std::sort(diagonal_segment.begin(),diagonal_segment.end(),DiagonalSeeds::cmp_subject_end);
    for (size_t i = 1; i < diagonal_segment.size(); i++)
    {
        int tmp = diagonal_segment[i].score; 
        int max_pre = i;
        for(j = 0; j < i; j++) 
        {
            int tdiff = static_cast<int>(diagonal_segment[j].j) - static_cast<int>(diagonal_segment[i].j) ;
            int qdiff = static_cast<int>(diagonal_segment[j].i) - static_cast<int>(diagonal_segment[i].i) ;
            int ref_gap = static_cast<int>(diagonal_segment[j].j) - static_cast<int>(diagonal_segment[i].j) - static_cast<int>(diagonal_segment[i].len);
            int qry_gap = static_cast<int>(diagonal_segment[j].i) - static_cast<int>(diagonal_segment[i].i) - static_cast<int>(diagonal_segment[i].len);

            if(tdiff <= 0 || qdiff <= 0 || qdiff > max_gap || tdiff > max_gap) 
            { // we often have multiple hits to the same target pos, so we can't let them chain with either other
                continue;
            }
            // if(ref_gap <= 0 || qry_gap <= 0 || qry_gap > max_gap || ref_gap > max_gap) 
            // { // we often have multiple hits to the same target pos, so we can't let them chain with either other
            //         continue;
            // }
            match_score = diagonal_segment[j].score;
            gap_cost = (diffdiff == 0 ? 0 : 0.01 * match_score * diffdiff + 0.5 * log2(diffdiff)); // affine gap penalty a la minimap2
            qdiff = (qdiff > tdiff ? tdiff : qdiff); // now not qdiff but the minimum difference
            score = diagonal_segment[j].score + (qdiff > match_score ? match_score : qdiff) - gap_cost; 
             if (tmp <= score)
             {  
                 tmp = score;
                 max_pre = j;
             } 
        }
        p1_score[i] = tmp;
        p1_track[i] = max_pre;
        if (tmp == diagonal_segment[i].score)
        {
            p1_track[i] = -1;
        }
    }
    
    bool *tracked = new bool [s_];
    memset(tracked, false, s_ * sizeof(bool));
    int best;
    bool flag = true;
    for (size_t j = s_ -1; j >= 0; j--)
    { 
        if (tracked[j])
        {
            continue;
        }
        int c = j;
        while (true)
        {
            tracked[c] = true;
            chained_seed.push_back(diagonal_segment[c]);
            best = p1_track[c];
            if (best == -1 || tracked[best])
            {
                break;
            }
            c = p1_track[c];
        }
    }

    delete [] tracked;
    delete [] p1_score;
    delete [] p1_track;

    return chained_seed;
}



#endif // __FINDSEEDSCHAIN_H__