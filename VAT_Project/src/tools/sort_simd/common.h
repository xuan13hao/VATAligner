#ifndef __COMMON_H__
#define __COMMON_H__



#include <iostream>
#include <x86intrin.h>
#include <string>
#include <cstdint>
#include <cassert>
#include <omp.h>

/**
 * Common definitions
 */

class Tuple
{
		public:
		Tuple():
			key (),
			value ()
		{ 

		}
		Tuple(uint64_t key, uint64_t value):
			key (key),
			value (value)
		{ 

		}
		bool operator<(const Tuple &rhs) const
		{ 
			return key < rhs.key; 
		}
		uint64_t	key;
		uint64_t		value;
} __attribute__((packed));////Tells the compiler to allocate the variable x at a 16-byte aligned memory address instead of the default 4-byte alignment.

#if __AVX2__
#define AVX2
#endif

#if __AVX512F__
#define AVX512
#endif

/**
 * Common utility functions
 */

/**
 * Templated function for memory aligned initialization
 * @tparam T: data type
 * @param ptr: reference to pointer needing init
 * @param N: size of data
 * @param alignment_size: alignment chunk size
 */
template <typename T>
void aligned_init(T* &ptr, size_t N, size_t alignment_size=64);

template <typename T>
void print_arr(T *arr, int i, int j, const std::string &tag="");

void print_kvarr(std::pair<int,int> *arr, int i, int j, const std::string &tag="");
void print_kvarr(int64_t *arr, int i, int j, const std::string &tag="");

#endif // __COMMON_H__