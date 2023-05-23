#include <iostream>
#include <vector>
#include <algorithm>
#include <assert.h>  
using namespace std;

struct DiagonalSeeds {
    int i, j, len, score;
    DiagonalSeeds(int query_pos, int subject_pos, int len, int score) :
        i(query_pos),
        j(subject_pos),
        len(len),
        score(score) {}
};

bool compare_anchors(const DiagonalSeeds& a1, const DiagonalSeeds& a2) {
    if (a1.score != a2.score) {
        return a1.score > a2.score;
    }
    if (a1.len != a2.len) {
        return a1.len > a2.len;
    }
    if (a1.i != a2.i) {
        return a1.i < a2.i;
    }
    return a1.j < a2.j;
}

void anchor_align(vector<DiagonalSeeds>& anchors, vector<vector<DiagonalSeeds>>& synteny_blocks) {
    // Sort anchors by score, length, i, j
    sort(anchors.begin(), anchors.end(), compare_anchors);

    // Loop through all anchors and group them into synteny blocks
    for (const auto& anchor : anchors) {
        bool found_synteny_block = false;
        for (auto& block : synteny_blocks) {
            // Check if anchor is compatible with the synteny block
            const DiagonalSeeds& last_anchor = block.back();
            if (anchor.i > last_anchor.i && anchor.j > last_anchor.j) {
                block.push_back(anchor);
                found_synteny_block = true;
                break;
            }
        }
        if (!found_synteny_block) {
            // Create a new synteny block with the current anchor
            synteny_blocks.push_back({anchor});
        }
    }

    // Order the synteny blocks
    sort(synteny_blocks.begin(), synteny_blocks.end(), [](const vector<DiagonalSeeds>& a, const vector<DiagonalSeeds>& b) {
        const DiagonalSeeds& a1 = a.front();
        const DiagonalSeeds& a2 = a.back();
        const DiagonalSeeds& b1 = b.front();
        const DiagonalSeeds& b2 = b.back();
        if (a1.i != b1.i) {
            return a1.i < b1.i;
        }
        if (a2.i != b2.i) {
            return a2.i < b2.i;
        }
        return a1.j < b1.j;
    });
}
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