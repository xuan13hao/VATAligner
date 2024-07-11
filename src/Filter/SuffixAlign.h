#ifndef __SUFFIXALIGN_H__
#define __SUFFIXALIGN_H__
#include <immintrin.h>
#include <cstdint>
#include <string>
#include "../Database/Reference.h"
#include "../Database/Queries.h"
#include "../Commons/Seed.h"
template<typename _val>
bool suffix_align( const _val *query , const _val* subject, std::string spaced_seed)
{
    if(spaced_seed != "null")
    {
        // Define the spaced seed pattern (example: 1111001111) 
        const size_t seed_length = spaced_seed.length();

        const _val* q = query + VATParameters::seed_len;
        const _val* s = subject + VATParameters::seed_len;

        size_t seed_index = 0;

        while (*q != AlphabetSet<_val>::PADDING_CHAR && *s != AlphabetSet<_val>::PADDING_CHAR && seed_index < seed_length)
        {
            if (spaced_seed[seed_index] == '1')
            {
                if (*q != *s)
                {
                    return false;
                }
            }

            // Move to the next position
            ++q;
            ++s;
            ++seed_index;
        }
        return true; // Return true if no mismatches found
    }
        else
    {
        return true;
    }

}



#endif // __SUFFIXALIGN_H__