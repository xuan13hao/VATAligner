#pragma once

#include "sort_util.h"
#include "merge_util.h"
#include "utils.h"
#include "../common.h"

#ifdef AVX512
namespace avx512{
  void SIMDSort(size_t N, int *&arr);
  void SIMDSort(size_t N, int64_t *&arr);
  void SIMDSort(size_t N, float *&arr);
  void SIMDSort(size_t N, double *&arr);
  void SIMDSort(size_t N, std::pair<int,int> *&arr);
  void SIMDOrderBy(std::pair<int, int> *&result_arr, size_t N, std::pair<int, int> *arr, int order_by=0);
  void SIMDSort(size_t N, std::pair<float, float> *&arr);
  void SIMDSort(size_t N, Tuple *&arr);
  void SIMDSort(size_t N, std::pair<int64_t ,int64_t> *&arr);
  void SIMDSort(size_t N, std::pair<double, double> *&arr);
};
#endif
