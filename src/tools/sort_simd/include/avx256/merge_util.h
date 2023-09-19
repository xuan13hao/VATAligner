#ifndef __MERGE_UTIL_H__
#define __MERGE_UTIL_H__


#include "Utils.h"
#include "../../common.h"

#ifdef AVX2

namespace avx2 {
  template <typename InType, typename RegType>
  void MergeRuns8(InType *&arr, size_t N);
  template <typename InType, typename RegType>
  void MaskedMergeRuns8(InType *&arr, size_t N);
  template <typename InType, typename RegType>
  void MergeRuns4(InType *&arr, size_t N);
  template <typename InType, typename RegType>
  void MaskedMergeRuns4(InType *&arr, size_t N);
  template <typename InType, typename RegType>
  void MergePass8(InType *&arr, InType *buffer, size_t N, unsigned int run_size);
  template <typename InType, typename RegType>
  void MaskedMergePass8(InType *&arr, InType *buffer, size_t N, unsigned int run_size);
  template <typename InType, typename RegType>
  void MergePass4(InType *&arr, InType *buffer, size_t N, unsigned int run_size);
  template <typename InType, typename RegType>
  void MaskedMergePass4(InType *&arr, InType *buffer, size_t N, unsigned int run_size);
};

#endif
#endif // __MERGE_UTIL_H__