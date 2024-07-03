#include <iostream>


std::vector<uint64_t> extractMinimizers(const DNA* sequence, size_t sequenceLength, size_t kmerSize, size_t windowSize) {
    std::vector<uint64_t> minimizers;
    uint64_t minSeed = std::numeric_limits<uint64_t>::max();

    for (size_t i = 0; i <= sequenceLength - kmerSize; i++) {
        uint64_t seed = 0;
        double f = 0;

        for (size_t j = 0; j < kmerSize; j++) {
            DNA l = sequence[i + j];
            if (l == AlphabetAttributes<DNA>::MASK_CHAR || l == AlphabetSet<DNA>::PADDING_CHAR) {
                seed = 0;
                break;
            }
            l = mask_critical(l);
            unsigned r = ReducedAlpha<DNA>::reduction(l);
            f += background_freq[r];
            seed *= ReducedAlpha<DNA>::reduction.size();
            seed += static_cast<uint64_t>(r);
        }

        seed = murmur_hash()(seed);
        if (use_seed_freq<DNA>() && f > VATParameters::max_seed_freq)
            seed = 0;

        if (seed != 0 && seed < minSeed)
            minSeed = seed;

        if (i % windowSize == windowSize - 1) {
            minimizers.push_back(minSeed);
            minSeed = std::numeric_limits<uint64_t>::max();
        }
    }

    return minimizers;
}

// Test function to print the minimizers
void printMinimizers(const std::vector<uint64_t>& minimizers) {
    std::cout << "Minimizers: ";
    for (const auto& minimizer : minimizers) {
        std::cout << minimizer << " ";
    }
    std::cout << std::endl;
}

int main() {
    // Test case
    const char* sequence = "ACGTACGTACGT";
    size_t sequenceLength = strlen(sequence);
    size_t kmerSize = 4;
    size_t windowSize = 4;

    // Call the extractMinimizers function
    std::vector<uint64_t> minimizers = extractMinimizers(sequence, sequenceLength, kmerSize, windowSize);

    // Print the resulting minimizers
    printMinimizers(minimizers);

    return 0;
}
