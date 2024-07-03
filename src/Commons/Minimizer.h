#ifndef __MINIMIZER_H__
#define __MINIMIZER_H__
/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param p      minimizers; p->a[i].x is the 2k-bit hash value;
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

vector<string> get_kmers(string sequence, int k) {
    vector<string> kmers;

    for (int i = 0; i <= sequence.length() - k; i++) {
        string kmer = sequence.substr(i, k);
        kmers.push_back(kmer);
    }

    return kmers;
}

string get_minimizer(string sequence, int k, int w) {
    vector<string> kmers = get_kmers(sequence, k);
    vector<string> window_kmers;
    string minimizer;

    for (int i = 0; i <= kmers.size() - w; i++) {
        window_kmers = vector<string>(kmers.begin() + i, kmers.begin() + i + w);
        sort(window_kmers.begin(), window_kmers.end());
        string window_minimizer = window_kmers[0];

        if (minimizer.empty() || window_minimizer < minimizer) {
            minimizer = window_minimizer;
        }
    }

    return minimizer;
}



vector<string> spaced_kmers(string sequence, string pattern, int k) {
    vector<string> kmers;
    int n = sequence.length();
    int num_positions = pattern.length();
    
    for (int i = 0; i <= n - k; i++) {
        bool valid_kmer = true;
        string kmer = "";
        for (int j = 0; j < num_positions; j++) {
            if (pattern[j] == '1') {
                if (sequence[i+j] == 'A' || sequence[i+j] == 'C' || sequence[i+j] == 'G' || sequence[i+j] == 'T') {
                    kmer += sequence[i+j];
                } else {
                    valid_kmer = false;
                    break;
                }
            } else {
                kmer += '-';
            }
        }
        if (valid_kmer) {
            kmers.push_back(kmer);
        }
    }
    return kmers;
}

int main() {
    string sequence = "AGCTACGTACGTA";
    string pattern = "1010";
    int k = 4;
    vector<string> kmers = spaced_kmers(sequence, pattern, k);
    for (string kmer : kmers) {
        cout << kmer << endl;
    }
    return 0;
}

#endif // __MINIMIZER_H__