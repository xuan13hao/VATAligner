#include <iostream>
#include <vector>

struct Seed {
    int x, y, l;
};

std::vector<Seed> diagonalSeedsChaining(const std::vector<Seed>& seeds, int maxGap) {
    std::vector<Seed> chains;
    std::vector<int> maxRight(seeds.size());
    int maxChainLen = 0;
    int maxChainId = -1;

    for (int i = seeds.size() - 1; i >= 0; --i) {
        maxRight[i] = i;
        for (int j = i + 1; j < seeds.size(); ++j) {
            int xGap = seeds[j].x - seeds[i].x;
            int yGap = seeds[j].y - seeds[i].y;
            if (xGap - yGap > maxGap) {
                break;
            }
            if (yGap - xGap > maxGap) {
                continue;
            }
            if (maxRight[j] > maxRight[i]) {
                maxRight[i] = maxRight[j];
            }
        }
        int chainLen = maxRight[i] - i + 1;
        if (chainLen >= maxChainLen) {
            maxChainLen = chainLen;
            maxChainId = i;
        }
    }

    if (maxChainId >= 0) {
        for (int i = maxChainId; i <= maxRight[maxChainId]; ++i) {
            chains.push_back(seeds[i]);
        }
    }

    return chains;
}

int main() {
    // Example usage
    std::vector<Seed> seeds = {{10, 20, 5}, {15, 25, 7}, {30, 40, 8}, {35, 45, 4}};
    int maxGap = 10;
    std::vector<Seed> chains = diagonalSeedsChaining(seeds, maxGap);

    std::cout << "Chains:\n";
    for (const auto& seed : chains) {
        std::cout << "x: " << seed.x << ", y: " << seed.y << ", l: " << seed.l << '\n';
    }

    return 0;
}
