#include <iostream>
#include <vector>
#include <algorithm>
#include "../basic/Seed.h"



template<typename _locr, typename _locl>
std::vector<DiagonalSeeds<_locr, _locl>> chainChimericSeeds(const std::vector<DiagonalSeeds<_locr, _locl>>& diagonalSeeds, int maxGap) {
    std::vector<DiagonalSeeds<_locr, _locl>> chimericSeeds;
    std::vector<DiagonalSeeds<_locr, _locl>> nonChimericSeeds;

    // Sort the diagonal seeds based on the query position
    std::vector<DiagonalSeeds<_locr, _locl>> sortedSeeds = diagonalSeeds;
    std::sort(sortedSeeds.begin(), sortedSeeds.end(),
              [](const DiagonalSeeds<_locr, _locl>& seed1, const DiagonalSeeds<_locr, _locl>& seed2) {
                  return seed1.i < seed2.i;
              });

    // Perform chaining algorithm
    std::vector<DiagonalSeeds<_locr, _locl>> currentChimera;
    int prevQueryPos = sortedSeeds[0].i;
    for (const auto& seed : sortedSeeds) {
        int queryPos = seed.i;
        int subjectPos = seed.j;

        // Check if the seed can be chained to the current chimera
        if (!currentChimera.empty() && queryPos - prevQueryPos <= maxGap) {
            // Check if the subject positions are different, indicating a split-read
            if (subjectPos != currentChimera.back().j) {
                currentChimera.push_back(seed);
            }
        } else {
            // Check if the current chimera is a valid chimeric seed
            if (currentChimera.size() > 1) {
                chimericSeeds.insert(chimericSeeds.end(), currentChimera.begin(), currentChimera.end());
            } else if (!currentChimera.empty()) {
                nonChimericSeeds.push_back(currentChimera.front());
            }
            currentChimera.clear();
            currentChimera.push_back(seed);
        }

        prevQueryPos = queryPos;
    }

    // Check if the last chimera is a valid chimeric seed
    if (currentChimera.size() > 1) {
        chimericSeeds.insert(chimericSeeds.end(), currentChimera.begin(), currentChimera.end());
    } else if (!currentChimera.empty()) {
        nonChimericSeeds.push_back(currentChimera.front());
    }

    // Merge the chimeric and non-chimeric seeds into a single vector
    std::vector<DiagonalSeeds<_locr, _locl>> allSeeds;
    allSeeds.reserve(chimericSeeds.size() + nonChimericSeeds.size());
    allSeeds.insert(allSeeds.end(), chimericSeeds.begin(), chimericSeeds.end());
    allSeeds.insert(allSeeds.end(), nonChimericSeeds.begin(), nonChimericSeeds.end());

    return allSeeds;
}

template<typename _locr, typename _locl>
std::vector<DiagonalSeeds<_locr, _locl>> findchimericSeeds(std::vector<DiagonalSeeds<_locr, _locl>>& seeds, int maxDistance) {
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

        for (int j = 0; j < i; ++j) {
            int distance = std::abs(static_cast<int>(seeds[i].j - seeds[j].j));
            // if (distance < maxDistance || !seeds[i].isChimericMapping(seeds[j]))
            //     continue;
            //weakly compatible
            if (!seeds[i].isChimericMapping(seeds[j]))
                continue;
            int score = dp[j] + seeds[i].score;
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
