#include <iostream>
#include <string>
#include <vector>

// Function to compute minimizers
std::vector<std::string> computeMinimizers(const std::string& sequence, int w, int k) {
    std::vector<std::string> minimizers;
    int sequenceLength = sequence.length();

    if (w > sequenceLength || k > w || w <= 0 || k <= 0) {
        // Invalid input parameters
        return minimizers;
    }

    for (int i = 0; i <= sequenceLength - w; i++) {
        std::string window = sequence.substr(i, w);

        // Initialize the smallest k-mer with the first k-mer in the window
        std::string smallestKmer = window.substr(0, k);

        for (int j = 1; j <= w - k; j++) {
            std::string currentKmer = window.substr(j, k);

            // Compare the current k-mer with the smallest k-mer
            if (currentKmer < smallestKmer) {
                smallestKmer = currentKmer;
            }
        }

        minimizers.push_back(smallestKmer);
    }

    return minimizers;
}

// Function to compute the hash value for a k-mer
std::size_t computeHashValue(const std::string& kmer) {
    std::hash<std::string> hashFunction;
    return hashFunction(kmer);
}


int main() {
    std::string sequence = "CCCCTTCGATGATTTCGCGCGTAAACACGTTGACGATTTTTCCGTTTTTGATGACGAGATCGGCCTTCGTTTGTTTTGCGGCGGCGGCGATTTGCTTCTGTAATGTCGAATGCATGTGGC"; // Your input sequence
    int w = 15; // Window size
    int k = 8;  // K-mer size

    std::vector<std::string> minimizers = computeMinimizers(sequence, w, k);

    std::cout << "Minimizers:" << std::endl;
    for (const std::string& minimizer : minimizers) {
        std::cout <<minimizer<<":"<< computeHashValue(minimizer) << std::endl;
    }

    return 0;
}
