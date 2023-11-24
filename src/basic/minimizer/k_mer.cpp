#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <set>
std::map<int, std::string> generateKmerMap(const std::string& sequence, int k) {
    std::map<int, std::string> kmerMap; // Key is the sum, and value is the k-mer
    int sequenceLength = sequence.length();

    // The total number of possible k-mers is 4^k
    int totalKmers = 1;
    for (int i = 0; i < k; ++i) {
        totalKmers *= 4;
    }

    for (int i = 0; i < totalKmers; ++i) {
        std::string kmer;

        // Generate a k-mer based on the current index 'i'
        int index = i;
        int sum = 0; // Initialize the sum for this k-mer

        for (int j = 0; j < k; ++j) {
            int nucleotideIndex = index % 4; // There are 4 nucleotides (A, C, G, T)
            char nucleotide;

            // Map nucleotides to numbers
            if (nucleotideIndex == 0) {
                nucleotide = 'A';
                sum = sum * 10 + 1;
            } else if (nucleotideIndex == 1) {
                nucleotide = 'C';
                sum = sum * 10 + 2;
            } else if (nucleotideIndex == 2) {
                nucleotide = 'G';
                sum = sum * 10 + 3;
            } else {
                nucleotide = 'T';
                sum = sum * 10 + 4;
            }

            kmer += nucleotide; // Build the k-mer
            index /= 4;
        }

        // Store the k-mer as a value, and the sum as the key
        kmerMap[sum] = kmer;
    }

    return kmerMap;
}

int computeKmerHash(const std::string& kmer) {
    int hash = 0;

    for (char nucleotide : kmer) {
        if (nucleotide == 'A') {
            hash = hash * 10 + 1;
        } else if (nucleotide == 'C') {
            hash = hash * 10 + 2;
        } else if (nucleotide == 'G') {
            hash = hash * 10 + 3;
        } else if (nucleotide == 'T') {
            hash = hash * 10 + 4;
        }
    }

    return hash;
}


// Function to represent a subsequence as a two-digit number
int subsequenceToKey(const std::string& subsequence) {
    // Define a map to map nucleotides to numbers
    std::map<char, int> nucleotideToNumber = {
        {'A', 1},
        {'C', 2},
        {'G', 3},
        {'T', 4}
    };

    int key = 0;
    for (char nucleotide : subsequence) {
        key = key * 10 + nucleotideToNumber[nucleotide];
    }
    return key;
}

// Function to cluster k-mers based on common subsequences
std::map<int, std::set<std::string>> clusterKmers(const std::map<int, std::string>& kmerMap, int k) {
    std::map<int, std::set<std::string>> clusters;

    for (const auto& entry : kmerMap) {
        const std::string& kmer = entry.second;

        // Iterate through the k-mer to find all subsequences of length (k-1)
        for (size_t i = 0; i <= kmer.length() - (k - 1); ++i) {
            std::string subsequence = kmer.substr(i, k - 1);
            int key = subsequenceToKey(subsequence);
            clusters[key].insert(kmer);
        }
    }

    return clusters;
}

int main() {
    std::string sequence = "ACGT"; // Replace with your DNA sequence
    int k = 3; // Size of the k-mer

    std::map<int, std::string> kmerMap = generateKmerMap(sequence, k);

    std::map<int, std::set<std::string>> clusters = clusterKmers(kmerMap, k);

    // Printing the clusters
    for (const auto& entry : clusters) {
        int key = entry.first;
        const std::set<std::string>& kmerSet = entry.second;

        std::cout << "Subsequence (Key): " << key << " K-mers: ";
        for (const std::string& kmer : kmerSet) {
            std::cout << kmer << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}