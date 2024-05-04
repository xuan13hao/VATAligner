#ifndef __CHAINING_H__
#define __CHAINING_H__

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

template<typename _locr, typename _locl>
std::vector<DiagonalSeeds<_locr, _locl>> chainingSeeds(std::vector<DiagonalSeeds<_locr, _locl>>& seeds, int maxDistance, int maxIndel, int maxLocalDistance) {
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
            int queryDistance = std::abs(static_cast<int>(seeds[i].i+ seeds[i].len- seeds[j].i - seeds[j].len));
            int targetDistance = std::abs(static_cast<int>(seeds[i].j + seeds[i].len - seeds[j].j - seeds[j].len));
            int indel = std::abs(static_cast<int>(queryDistance - targetDistance));


            bool colinear = (seeds[i].i + seeds[i].len < seeds[j].i + seeds[j].len && seeds[i].j + seeds[i].len < seeds[j].j + seeds[j].len) || (seeds[i].i + seeds[i].len> seeds[j].i + seeds[j].len && seeds[i].j + seeds[i].len> seeds[j].j+ seeds[j].len);
            bool nonOverlapping = (seeds[i].i + seeds[i].len> seeds[j].i + seeds[j].len) || (seeds[i].j  + seeds[i].len> seeds[j].j + seeds[j].len);
            bool localDistance = (queryDistance < maxLocalDistance && targetDistance < maxLocalDistance);
            bool smallIndel = (indel < maxIndel);
            //strongly and weakly compatible
            // if (queryDistance > maxDistance || targetDistance > maxDistance || !colinear || !nonOverlapping || !localDistance || !smallIndel)
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


#endif // __CHAINING_H__