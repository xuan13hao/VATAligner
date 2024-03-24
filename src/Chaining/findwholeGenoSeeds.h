#ifndef __FINDWHOLEGENOSEEDS_H__
#define __FINDWHOLEGENOSEEDS_H__

#include <vector>
#include <algorithm> 
#include <string>
#include <cmath>
#include <assert.h>
#include "../basic/Seed.h"
// #include "kvec.h"
using std::vector;
/*
template<typename _locr, typename _locl>
vector<DiagonalSeeds<_locr,_locl> > findWholeGenSeeds(vector<DiagonalSeeds<_locr,_locl> > & diagonal_segment,int max_gap)
{
    vector<DiagonalSeeds<_locr,_locl> > chained_seed;
    int s_ =  diagonal_segment.size();
    vector<int> dp(s_+1, 0);
    vector<int> pre(s_+1, -1);
    dp[0] = diagonal_segment[0].score;
    std::sort(diagonal_segment.begin(),diagonal_segment.end(),DiagonalSeeds<_locr,_locl>::cmp_subject_end);
    
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
    vector<vector<DiagonalSeeds<_locr, _locl>>> synteny_blocks;
    anchor_align(chained_seed, synteny_blocks);
    // Concatenate synteny blocks into results vector
    vector<DiagonalSeeds<_locr, _locl>> results;
    for (const auto& block : synteny_blocks) {
        results.insert(results.end(), block.begin(), block.end());
    }
    return results;

}
*/


template<typename _locr, typename _locl>
std::vector<DiagonalSeeds<_locr, _locl>> ChainWGSSeeds(std::vector<DiagonalSeeds<_locr, _locl>>& seeds, int maxDistance, int maxIndel, int maxLocalDistance) {
    std::vector<int> dp(seeds.size());
    std::vector<int> prev(seeds.size(), -1);
    std::vector<int> maxLen(seeds.size());
    std::vector<int> maxIdx(seeds.size());

    int bestScore = 0;
    int bestIdx = -1;

    for (int i = 0; i < seeds.size(); ++i) {
        dp[i] = seeds[i].score;
        maxLen[i] = seeds[i].len;
        maxIdx[i] = i;
//i - j <= 50
        for (int j = std::max(0, i - 50); j < i; ++j) {
        // for (int j = 0; j < i; ++j) {
            int queryDistance = std::abs(static_cast<int>(seeds[i].i+ seeds[i].len- seeds[j].i - seeds[j].len));
            int targetDistance = std::abs(static_cast<int>(seeds[i].j + seeds[i].len - seeds[j].j - seeds[j].len));
            int indel = std::abs(static_cast<int>(queryDistance - targetDistance));

            bool localDistance = (queryDistance < maxLocalDistance && targetDistance < maxLocalDistance);
            bool smallIndel = (indel < maxIndel);
            // weakly compatible
            if (queryDistance > maxDistance || targetDistance > maxDistance)
                continue;
            // if (queryDistance > maxDistance || targetDistance > maxDistance || !localDistance || !smallIndel)
            //     continue;

            int score = dp[j] + seeds[i].len;
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

    // vector<vector<DiagonalSeeds<_locr, _locl>>> synteny_blocks;
    // anchor_align(chainedSeeds, synteny_blocks);
    // // Concatenate synteny blocks into results vector
    // vector<DiagonalSeeds<_locr, _locl>> results;
    // for (const auto& block : synteny_blocks) {
    //     results.insert(results.end(), block.begin(), block.end());
    // }
    // return results;

template<typename _locr, typename _locl>
std::vector<DiagonalSeeds<_locr, _locl>> findWholeGenSeeds(std::vector<DiagonalSeeds<_locr, _locl>>& seeds, int maxGap)
{
    std::vector<int> dp(seeds.size());
    std::vector<int> prev(seeds.size(), -1);
    std::vector<int> maxLen(seeds.size());
    std::vector<int> maxIdx(seeds.size());

    int bestScore = 0;
    int bestIdx = -1;
    for (int i = 0; i < seeds.size(); ++i) {
        dp[i] = seeds[i].len;
        maxLen[i] = seeds[i].len;
        maxIdx[i] = i;

        for (int j = 0; j < i; ++j) {
            int gap = seeds[i].i + seeds[i].len- seeds[j].i - seeds[j].len;
            if (gap > maxGap)
                continue;

            int score = dp[j] + seeds[i].len;
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

    vector<vector<DiagonalSeeds<_locr, _locl>>> synteny_blocks;
    anchor_align(chainedSeeds, synteny_blocks);
    // Concatenate synteny blocks into results vector
    vector<DiagonalSeeds<_locr, _locl>> results;
    for (const auto& block : synteny_blocks) {
        results.insert(results.end(), block.begin(), block.end());
    }
    return results;
}

template<typename _locr, typename _locl>
void anchor_align(vector<DiagonalSeeds<_locr, _locl>>& anchors, vector<vector<DiagonalSeeds<_locr, _locl>>>& synteny_blocks) {
    // Sort anchors by score, length, i, j
    sort(anchors.begin(), anchors.end(), [](const DiagonalSeeds<_locr, _locl>& a, const DiagonalSeeds<_locr, _locl>& b) {
        if (a.score != b.score) {
            return a.score > b.score;
        }
        if (a.len != b.len) {
            return a.len > b.len;
        }
        if (a.i != b.i) {
            return a.i < b.i;
        }
        return a.j < b.j;
    });

    // Loop through all anchors and group them into synteny blocks
    for (const auto& anchor : anchors) {
        bool found_synteny_block = false;
        for (auto& block : synteny_blocks) {
            // Check if anchor is compatible with the synteny block
            const DiagonalSeeds<_locr, _locl>& last_anchor = block.back();
            if (anchor.i > last_anchor.i && anchor.j > last_anchor.j) {
                block.push_back(anchor);
                found_synteny_block = true;
                break;
            }
        }
        if (!found_synteny_block) {
            // Create a new synteny block with the current anchor
            synteny_blocks.push_back({ anchor });
        }
    }

    // Order the synteny blocks
    sort(synteny_blocks.begin(), synteny_blocks.end(), [](const vector<DiagonalSeeds<_locr, _locl>>& a, const vector<DiagonalSeeds<_locr, _locl>>& b) {
        const DiagonalSeeds<_locr, _locl>& a1 = a.front();
        const DiagonalSeeds<_locr, _locl>& a2 = a.back();
        const DiagonalSeeds<_locr, _locl>& b1 = b.front();
        const DiagonalSeeds<_locr, _locl>& b2 = b.back();
        if (a1.i != b1.i) {
            return a1.i < b1.i;
        }
        if (a2.i != b2.i) {
            return a2.i < b2.i;
        }
        return a1.j < b1.j;
    });
}

/*
int main() {
    vector<DiagonalSeeds> anchors = {
        DiagonalSeeds(1, 3, 10, 50),
        DiagonalSeeds(2, 5, 8, 45),
        DiagonalSeeds(3, 7, 6, 40),
        DiagonalSeeds(4, 1, 7, 30),
        DiagonalSeeds(5, 2, 9, 35),
        DiagonalSeeds(6, 4, 5, 40),
        DiagonalSeeds(7, 6, 12, 60),
        DiagonalSeeds(1, 5, 11, 61),
        DiagonalSeeds(2, 6, 14, 48)
    };

    vector<vector<DiagonalSeeds>> synteny_blocks;
    anchor_align(anchors, synteny_blocks);

    for (const auto& block : synteny_blocks) {
        cout << "Synteny block:" << endl;
        for (const auto& anchor : block) {
            cout << "i: " << anchor.i << " j: " << anchor.j << " len: " << anchor.len << " score: " << anchor.score << endl;
        }
    }

    return 0;
}
*/

#endif // __FINDWHOLEGENOSEEDS_H__