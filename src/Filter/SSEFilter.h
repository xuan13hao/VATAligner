#ifndef __SSEFILTER_H__
#define __SSEFILTER_H__



#include <immintrin.h>  // For AVX2 and SIMD instructions
#include <iostream>
#ifdef __SSSE3__
#include <tmmintrin.h>
#endif
#include <stdint.h>
#ifdef _MSC_VER
#include <intrin.h>
#endif
#include "../Commons/ReducedAlpha.h"

unsigned popcount64(unsigned long long x)
{
#ifdef _MSC_VER
    // Use MSVC's intrinsic function to count bits
    return (unsigned)__popcnt64(x);
#else
    // Use GCC/Clang's built-in function to count bits
    return __builtin_popcountll(x);
#endif
}

//The number of 1 bits in the value of x
unsigned popcount_3(uint64_t x)
{
	const uint64_t m1  = 0x5555555555555555; //binary: 0101...
	const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
	const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
	const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits
    return (x * h01)>>56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...
}

template<typename _val>
unsigned match_block(const _val *x, const _val *y)
{
	static const __m128i mask = _mm_set1_epi8(0x7F);
	__m128i r1 = _mm_loadu_si128 ((__m128i const*)(x));
	__m128i r2 = _mm_loadu_si128 ((__m128i const*)(y));
	r2 = _mm_and_si128(r2, mask);
	return _mm_movemask_epi8(_mm_cmpeq_epi8(r1, r2));
}
template<typename _val>
unsigned match_block_avx2(const _val *x, const _val *y)
{
    static const __m256i mask = _mm256_set1_epi8(0x7F);
    __m256i r1 = _mm256_loadu_si256((__m256i const*)(x));
    __m256i r2 = _mm256_loadu_si256((__m256i const*)(y));
    r2 = _mm256_and_si256(r2, mask);

    // Perform comparison using 256-bit AVX2 registers
    __m256i cmp_result = _mm256_cmpeq_epi8(r1, r2);

    // Extract a bitmask indicating which bytes are equal
    // int bitmask = _mm256_movemask_epi8(cmp_result);

    return _mm256_movemask_epi8(cmp_result);
}


// Function to count the number of set bits (popcount) in an integer
unsigned popcount3(unsigned x) {
    unsigned count = 0;
    while (x) {
        count += x & 1;
        x >>= 1;
    }
    return count;
}

// Function to count the number of set bits in a 256-bit AVX2 register
int popcount_256(__m256i vec) {
    int total_popcount = 0;

    // Store the 256-bit vector into an array of 32 bytes
    alignas(32) unsigned char bytes[32];
    _mm256_storeu_si256((__m256i*)bytes, vec);

    // Count the set bits in each byte
    for (int i = 0; i < 32; ++i) {
        total_popcount += popcount3(bytes[i]);
    }

    return total_popcount;
}

// Function to compute the Hamming distance between two sequences using AVX2
template<typename _val>
int hamming_distance_avx2(const _val *q, const _val *s, size_t length) {
    size_t i = 0;
    int hamming_dist = 0;

    // Process 32 bytes at a time using AVX2
    for (; i + 31 < length; i += 32) {
        __m256i vec_q = _mm256_loadu_si256((__m256i const*)(q + i));  // Load 32 bytes from q
        __m256i vec_s = _mm256_loadu_si256((__m256i const*)(s + i));  // Load 32 bytes from s

        __m256i xor_result = _mm256_xor_si256(vec_q, vec_s);  // XOR the two blocks
        hamming_dist += popcount_256(xor_result);  // Count the number of differing bits
    }
	// AlphabetAttributes<_val>::ALPHABET[q[i]]
    // Process any remaining bytes
    for (; i < length; ++i) {
        hamming_dist += popcount_3(AlphabetAttributes<_val>::ALPHABET[q[i]] ^ AlphabetAttributes<_val>::ALPHABET[s[i]]);
    }

    return hamming_dist;
}
template<typename _val>
unsigned fast_match(const _val *q, const _val *s)
{ 
	return popcount_3(match_block(q-8, s-8)<<16 | match_block(q+8, s+8)); 
}

template<typename _val>
unsigned fast_match_short(const _val *q, const _val *s)
{
	// unsigned count = popcount_3((match_block_avx2(q-8, s-8) << 16) | match_block_avx2(q+8, s+8));
	// cout<<"count = "<<count<<"\n";
    return popcount_3((match_block_avx2(q-8, s-8) << 16) | match_block_avx2(q+8, s+8));
}

template<typename _val>
unsigned fast_match_forward(const _val *q, const _val *s)
{
    // return popcount64(match_block_avx2(q+16, s+16));
	return popcount_3(match_block_avx2(q+16, s+16));
}

template<typename _val>
unsigned fast_match_long(const _val *q, const _val *s)
{
    return popcount_3(((uint64_t)match_block_avx2(q-16, s-16) << 32) | match_block_avx2(q+16, s+16));
}

template<typename _val>
__m128i reduce_seq_ssse3(const __m128i &seq)
{
#ifdef __SSSE3__
	const __m128i *row = reinterpret_cast<const __m128i*>(ReducedAlpha<_val>::reduction.map8());
	__m128i high_mask = _mm_slli_epi16(_mm_and_si128(seq, _mm_set1_epi8(0x10)), 3);
	__m128i seq_low = _mm_or_si128(seq, high_mask);
	__m128i seq_high = _mm_or_si128(seq, _mm_xor_si128(high_mask, _mm_set1_epi8(0x80)));

	__m128i r1 = _mm_load_si128(row);
	__m128i r2 = _mm_load_si128(row+1);
	__m128i s1 = _mm_shuffle_epi8(r1, seq_low);
	__m128i s2 = _mm_shuffle_epi8(r2, seq_high);
	return _mm_or_si128(s1, s2);
#endif
}

template<typename _val>
__m128i reduce_seq_generic(const __m128i &seq)
{
	__m128i r;
	_val* s = (_val*)&seq;
	uint8_t* d = (uint8_t*)&r;
	for(unsigned i=0;i<16;++i)
		*(d++) = ReducedAlpha<_val>::reduction(*(s++));
	return r;
}

template<typename _val>
__m128i reduce_seq(const __m128i &seq)
{
	if(VATParameters::have_ssse3) {
#ifdef __SSSE3__
		return reduce_seq_ssse3<_val>(seq);
#else
		return reduce_seq_generic<_val>(seq);
#endif
	} else
		return reduce_seq_generic<_val>(seq);
}

template<typename _val>
unsigned match_block_reduced(const _val *x, const _val *y)
{
	static const __m128i mask = _mm_set1_epi8(0x7F);
	__m128i r1 = _mm_loadu_si128 ((__m128i const*)(x));
	__m128i r2 = _mm_loadu_si128 ((__m128i const*)(y));
	r2 = _mm_and_si128(r2, mask);
	r1 = reduce_seq<_val>(r1);
	r2 = reduce_seq<_val>(r2);
	return _mm_movemask_epi8(_mm_cmpeq_epi8(r1, r2));
}

template<typename _val>
uint64_t reduced_match32(const _val* q, const _val *s, unsigned len)
{
	uint64_t x = match_block_reduced(q+16, s+16)<<16 | match_block_reduced(q,s);
	if(len < 32)
		x &= (1 << len) - 1;
	return x;
}


#endif // __SSEFILTER_H__