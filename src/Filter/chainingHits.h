#ifndef __CHAININGHITS_H__
#define __CHAININGHITS_H__

#include <iostream>
#include <vector>
#include <algorithm>

template<typename _locr, typename _locl>
vector<vector<Hits<_locr, _locl>>> chaining_algorithm(const vector<Hits<_locr, _locl>>& hits, int max_gap)
{
    vector<vector<Hits<_locr, _locl>>> chains;
    vector<bool> used(hits.size(), false); // Track if a hit is already used in a chain

    // Sort hits based on subject values
    vector<Hits<_locr, _locl>> sortedHits = hits;
    sort(sortedHits.begin(), sortedHits.end(), Hits<_locr, _locl>::cmp_subject);

    // Iterate through the sorted hits
    for (int i = 0; i < sortedHits.size(); i++)
    {
        if (used[i])
            continue;

        vector<Hits<_locr, _locl>> chain;
        chain.push_back(sortedHits[i]);
        used[i] = true;

        // Extend the chain in both directions
        for (int j = i + 1; j < sortedHits.size(); j++)
        {
            if (used[j])
                continue;

            const Hits<_locr, _locl>& prevHit = chain.back();
            const Hits<_locr, _locl>& currHit = sortedHits[j];

            // Check if the hits can be chained based on the max_gap value
            if (currHit.subject_ - prevHit.subject_ <= max_gap)
            {
                chain.push_back(currHit);
                used[j] = true;
            }
            else
            {
                break; // Hits are too far apart, stop extending the chain
            }
        }

        chains.push_back(chain);
    }

    return chains;
}


int main()
{
    // Example usage
    std::vector<Hits<unsigned, int>> hits;
    // Populate the hits vector with actual hits
    
    // Call the chaining algorithm
    chaining_algorithm(hits);
    
    return 0;
}
#endif // __CHAININGHITS_H__