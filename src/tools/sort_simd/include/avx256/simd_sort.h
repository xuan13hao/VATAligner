#ifndef __SIMD_SORT_H__
#define __SIMD_SORT_H__


#include "sort_util.h"
#include "merge_util.h"
#include "Utils.h"
#include "../../common.h"

#ifdef AVX2
namespace avx2{
  void SIMDSort(size_t N, Tuple *&arr);
  void SIMDSort(size_t N, std::pair<uint64_t ,uint64_t> *&arr);
};
#endif
#endif // __SIMD_SORT_H__