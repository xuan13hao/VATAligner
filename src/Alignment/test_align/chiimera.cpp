#include <iostream>
#include <vector>
#include <algorithm>
// Define the DiagonalSeeds class template
template<typename _locr, typename _locl>
class DiagonalSeeds {
public:
    _locr i;         // Query position
    _locl j;         // Subject position
    int len;        // Length
    int score;      // Score

    DiagonalSeeds(_locr query_pos, _locl subject_pos, int length, int score)
        : i(query_pos),
          j(subject_pos),
          len(length),
          score(score)
    {}
};


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


int main() {
    // Create a vector of DiagonalSeeds
    std::vector<DiagonalSeeds<int, int>> diagonalSeeds;
    diagonalSeeds.emplace_back(DiagonalSeeds<int, int>(10, 20, 5, 50));
    diagonalSeeds.emplace_back(DiagonalSeeds<int, int>(15, 25, 4, 40));
    diagonalSeeds.emplace_back(DiagonalSeeds<int, int>(30, 40, 7, 70));
    diagonalSeeds.emplace_back(DiagonalSeeds<int, int>(40, 50, 6, 60));
    diagonalSeeds.emplace_back(DiagonalSeeds<int, int>(55, 65, 5, 50));
    diagonalSeeds.emplace_back(DiagonalSeeds<int, int>(60, 70, 8, 80));

    // Specify the maximum gap for chaining
    int maxGap = 10;

    // Call the chainChimericSeeds function
    std::vector<DiagonalSeeds<int, int>> result = chainChimericSeeds(diagonalSeeds, maxGap);

    // Output the chimeric seeds
    std::cout << "Chimeric Seeds:" << std::endl;
    for (const auto& seed : result) {
        std::cout << "Query Pos: " << seed.i << ", Subject Pos: " << seed.j
                  << ", Length: " << seed.len << ", Score: " << seed.score << std::endl;
    }

    return 0;
}