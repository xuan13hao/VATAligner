#include <iostream>
#include <vector>
#include <algorithm>

template<typename _locr, typename _locl>
class DiagonalSeeds
{
public:
    DiagonalSeeds() :
        len(0)
    {}

    DiagonalSeeds(_locr query_pos, _locl subject_pos, int len, int score) :
        i(query_pos),
        j(subject_pos),
        len(len),
        score(score)
    {}

    // Other member functions or variables can be added here

    _locr i;      // Type of query position
    _locl j;      // Type of subject position
    int len;      // Length of the diagonal seed
    int score;    // Score associated with the diagonal seed
};

template<typename _locr, typename _locl>
std::vector<DiagonalSeeds<_locr, _locl>> chainSeeds(const std::vector<DiagonalSeeds<_locr, _locl>>& seeds, int maxGap, int maxChain)
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
            int gap = seeds[i].i - seeds[j].i - seeds[j].len;
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
    return chainedSeeds;
}

int main()
{
    // Example usage of chainSeeds with DiagonalSeeds

    std::vector<DiagonalSeeds<int, int>> seeds = {
        {10, 20, 5, 100},
        {15, 25, 10, 200},
        {30, 40, 8, 150},
        {45, 55, 6, 180},
        {60, 70, 7, 220}
    };

    int maxGap = 150;
    int maxChain = 3;

    std::vector<DiagonalSeeds<int, int>> chainedSeeds = chainSeeds(seeds, maxGap, maxChain);

    std::cout << "Chained Seeds:\n";
    for (const auto& seed : chainedSeeds) {
        std::cout << "i: " << seed.i << ", j: " << seed.j << ", len: " << seed.len << ", score: " << seed.score << "\n";
    }

    return 0;
}
