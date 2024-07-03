
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

template<typename _locr, typename _locl>
class DiagonalSeeds {
public:
    DiagonalSeeds() :
        len(0)
    {}

    DiagonalSeeds(_locr query_pos, _locl subject_pos, _locr len, _locl score) :
        i(query_pos),
        j(subject_pos),
        len(len),
        score(score)
    {}

    _locr i; // Query position
    _locl j; // Subject position
    _locr len; // Length
    _locl score; // Score
};

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
};

template<typename _locr, typename _locl>
void test_anchor_align() {
    // Create some test anchors
    vector<DiagonalSeeds<_locr, _locl>> anchors;
    anchors.push_back(DiagonalSeeds<_locr, _locl>(1, 2, 10, 100));
    anchors.push_back(DiagonalSeeds<_locr, _locl>(3, 4, 8, 90));
    anchors.push_back(DiagonalSeeds<_locr, _locl>(2, 3, 12, 95));
    anchors.push_back(DiagonalSeeds<_locr, _locl>(4, 5, 6, 80));
    anchors.push_back(DiagonalSeeds<_locr, _locl>(6, 7, 9, 85));

    // Perform anchor alignment
    vector<vector<DiagonalSeeds<_locr, _locl>>> synteny_blocks;
    anchor_align(anchors, synteny_blocks);

    // Print the resulting synteny blocks
    for (const auto& block : synteny_blocks) {
        for (const auto& anchor : block) {
            cout << "i: " << anchor.i << ", j: " << anchor.j << ", len: " << anchor.len << ", score: " << anchor.score << endl;
        }
        cout << "---------------------------" << endl;
    }
}

int main() {
    test_anchor_align<int, int>();  // Test with int as the template arguments for _locr and _locl

    return 0;
}

