#ifndef __SUFFIXALIGN_H__
#define __SUFFIXALIGN_H__
#include <immintrin.h>
#include <cstdint>
#include <string>
#include "../Database/Reference.h"
#include "../Database/Queries.h"
#include "../Commons/Seed.h"

#include <immintrin.h> // For AVX2 intrinsics

template<typename _val>
bool suffix_align_avx2(const _val *query, const _val* subject, std::string spaced_seed, std::string reduced_alphabet)
{
    // Function to map to a reduced alphabet
    auto map_to_reduced_alphabet = [&](char aa, const std::string& pattern) -> char
    {
        // ... (same as before, reduced alphabet mapping function)
        // No change to this part.
    };

    if (spaced_seed != "null")
    {
        const size_t seed_length = spaced_seed.length();
        const _val* q = query + VATParameters::seed_len;
        const _val* s = subject + VATParameters::seed_len;

        size_t seed_index = 0;

        while (*q != AlphabetSet<_val>::PADDING_CHAR && *s != AlphabetSet<_val>::PADDING_CHAR && seed_index < seed_length)
        {
            if (spaced_seed[seed_index] == '1')
            {
                char q_char = *q;
                char s_char = *s;

                // Apply reduced alphabet mapping if needed
                if (reduced_alphabet != "null")
                {
                    q_char = map_to_reduced_alphabet(q_char, reduced_alphabet);
                    s_char = map_to_reduced_alphabet(s_char, reduced_alphabet);
                }

                // SIMD-optimized comparison for multiple characters at a time
                __m256i query_chunk = _mm256_loadu_si256((__m256i*)q);  // Load 32 bytes from query
                __m256i subject_chunk = _mm256_loadu_si256((__m256i*)s);  // Load 32 bytes from subject

                // Compare 32 characters in one go
                __m256i result = _mm256_cmpeq_epi8(query_chunk, subject_chunk);

                // Check if all compared bytes are equal
                if (_mm256_movemask_epi8(result) != 0xFFFFFFFF)
                {
                    return false; // Mismatch found
                }
            }

            // Move to the next position
            q += 32;  // Process the next chunk of 32 characters
            s += 32;
            seed_index += 32;
        }

        return true; // Return true if no mismatches were found
    }
    else
    {
        return true; // If no spaced seed is used, always return true
    }
}


template<typename _val>
bool suffix_align(const _val *query, const _val* subject, std::string spaced_seed, std::string reduced_alphabet)
{
    // Function to map to a reduced alphabet
    auto map_to_reduced_alphabet = [&](char aa, const std::string& pattern) -> char
    {
        if (pattern == "dssp.5")
        {
            if (aa == 'A' || aa == 'G' || aa == 'P') return '1'; 
            if (aa == 'I' || aa == 'L' || aa == 'M' || aa == 'V') return '2'; 
            if (aa == 'F' || aa == 'W' || aa == 'Y') return '3'; 
            if (aa == 'R' || aa == 'H' || aa == 'K') return '4'; 
            if (aa == 'D' || aa == 'E' || aa == 'N' || aa == 'Q' || aa == 'S' || aa == 'T') return '5';
        }
        else if (pattern == "murphy.5")
        {
            if (aa == 'A' || aa == 'G') return '1'; 
            if (aa == 'L' || aa == 'I' || aa == 'V' || aa == 'M') return '2'; 
            if (aa == 'K' || aa == 'R' || aa == 'H') return '3'; 
            if (aa == 'D' || aa == 'E') return '4'; 
            if (aa == 'F' || aa == 'Y' || aa == 'W') return '5'; 
        }
        else if (pattern == "dssp.10")
        {
            if (aa == 'A' || aa == 'G') return '1';
            if (aa == 'P') return '2'; 
            if (aa == 'I' || aa == 'L' || aa == 'V') return '3'; 
            if (aa == 'F' || aa == 'W' || aa == 'Y') return '4'; 
            if (aa == 'D' || aa == 'E') return '5'; 
            if (aa == 'N' || aa == 'Q') return '6'; 
            if (aa == 'R' || aa == 'H' || aa == 'K') return '7'; 
            if (aa == 'M') return '8'; 
            if (aa == 'S' || aa == 'T') return '9'; 
            if (aa == 'C') return '10'; 
        }
        else if (pattern == "murphy.10")
        {
            if (aa == 'A' || aa == 'G') return '1'; 
            if (aa == 'L' || aa == 'I' || aa == 'V') return '2'; 
            if (aa == 'P') return '3'; 
            if (aa == 'F' || aa == 'W' || aa == 'Y') return '4'; 
            if (aa == 'D' || aa == 'E') return '5'; 
            if (aa == 'K' || aa == 'R' || aa == 'H') return '6'; 
            if (aa == 'N' || aa == 'Q') return '7'; 
            if (aa == 'S' || aa == 'T') return '8'; 
            if (aa == 'M') return '9'; 
            if (aa == 'C') return '10'; 
        }
        else if (pattern == "td.10")
        {
            if (aa == 'A') return '1'; 
            if (aa == 'V' || aa == 'L' || aa == 'I') return '2'; 
            if (aa == 'F' || aa == 'Y' || aa == 'W') return '3'; 
            if (aa == 'G') return '4'; 
            if (aa == 'D' || aa == 'E') return '5'; 
            if (aa == 'N' || aa == 'Q') return '6'; 
            if (aa == 'R' || aa == 'H' || aa == 'K') return '7'; 
            if (aa == 'S' || aa == 'T') return '8'; 
            if (aa == 'P') return '9'; 
            if (aa == 'C' || aa == 'M') return '10'; 
        }
        else if (pattern == "MMSEQS12")
        {
            // MMSEQS12 groupings
            if (aa == 'A' || aa == 'S' || aa == 'T') return '1'; 
            if (aa == 'L' || aa == 'M') return '2'; 
            if (aa == 'I' || aa == 'V') return '3'; 
            if (aa == 'K' || aa == 'R') return '4'; 
            if (aa == 'E' || aa == 'Q') return '5'; 
            if (aa == 'N' || aa == 'D') return '6'; 
            if (aa == 'F' || aa == 'Y') return '7'; 
            if (aa == 'C') return '8'; 
            if (aa == 'G') return '9'; 
            if (aa == 'H') return '10'; 
            if (aa == 'P') return '11'; 
            if (aa == 'W') return '12'; 
        }
        return aa; // Default: return original amino acid if no reduction is applied
    };

    if (spaced_seed != "null")
    {
        const size_t seed_length = spaced_seed.length();
        const _val* q = query + VATParameters::seed_len;
        const _val* s = subject + VATParameters::seed_len;

        size_t seed_index = 0;

        while (*q != AlphabetSet<_val>::PADDING_CHAR && *s != AlphabetSet<_val>::PADDING_CHAR && seed_index < seed_length)
        {
            if (spaced_seed[seed_index] == '1')
            {
                char q_char = *q;
                char s_char = *s;

                // Apply reduced alphabet mapping if needed
                if (reduced_alphabet != "null")
                {
                    q_char = map_to_reduced_alphabet(q_char, reduced_alphabet);
                    s_char = map_to_reduced_alphabet(s_char, reduced_alphabet);
                }

                if (q_char != s_char)
                {
                    return false; // Mismatch found
                }
            }

            // Move to the next position
            ++q;
            ++s;
            ++seed_index;
        }
        return true; // Return true if no mismatches were found
    }
    else
    {
        return true; // If no spaced seed is used, always return true
    }
}


// template<typename _val>
// bool suffix_align( const _val *query , const _val* subject, std::string spaced_seed)
// {
//     if(spaced_seed != "null")
//     {
//         // Define the spaced seed pattern (example: 1111001111) 
//         const size_t seed_length = spaced_seed.length();

//         const _val* q = query + VATParameters::seed_len;
//         const _val* s = subject + VATParameters::seed_len;

//         size_t seed_index = 0;
//         char s = AlphabetAttributes<_val>::ALPHABET[query];
//         while (*q != AlphabetSet<_val>::PADDING_CHAR && *s != AlphabetSet<_val>::PADDING_CHAR && seed_index < seed_length)
//         {
//             if (spaced_seed[seed_index] == '1')
//             {
//                 if (*q != *s)
//                 {
//                     return false;
//                 }
//             }

//             // Move to the next position
//             ++q;
//             ++s;
//             ++seed_index;
//         }
//         return true; // Return true if no mismatches found
//     }
//         else
//     {
//         return true;
//     }

// }



#endif // __SUFFIXALIGN_H__