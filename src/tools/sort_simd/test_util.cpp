#include "./avx256/utils.h"
#include "./avx256/simd_sort.h"
#include "common.h"
#include <random>
#define LO -100
#define HI 100
#define NNUM 4
//g++ test.cpp common.cpp avx256/* -march=native -o a
/**
 template<typename InType, typename RegType>
void LoadReg(RegType &r, InType *arr) {
  r = *((RegType *) arr);
}

template void LoadReg<uint64_t, __m256i>(__m256i &r, uint64_t *arr);
*/

  template <typename T>
  static void RandGenInt(T* &arr, size_t N, T lo, T hi) 
  {
    aligned_init<T>(arr, N);
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<T> dis(lo, hi);
    for(size_t i = 0; i < N; i++) 
    {
      arr[i] = dis(gen);
    //   cout<<"arr = "<<arr[i]<<endl;
    }
  }
// Load/Stores
template <typename InType, typename RegType>
void LoadReg(RegType &r, InType* arr);
template <typename InType, typename RegType>
void StoreReg(const RegType &r, InType* arr);
template<typename InType, typename RegType>
void StoreReg(const RegType &r, InType *arr) {
  *((RegType *) arr) = r;
}
template<typename InType, typename RegType>
void LoadReg(RegType &r, InType *arr) {
  r = *((RegType *) arr);
}


  int main(int argc, const char** argv) 
  {

    uint64_t *a;
    uint64_t *b;
    aligned_init<uint64_t>(a, 4);
    aligned_init<uint64_t>(b, 4);
    RandGenInt<uint64_t>(a, 4, -10, 10);
    RandGenInt<uint64_t>(b, 4, -10, 10);
    __m256i ra, rb;
    LoadReg(ra, a);
    LoadReg(rb, b);
    // MinMax4(ra, rb);
    StoreReg(ra, a);
    StoreReg(rb, b);
      // RandGenInt(rand_arr, N, lo, hi);
      // // aligned_init<int>(soln_arr, N);
      // std::copy(rand_arr, rand_arr + N, soln_arr);
      // for (size_t i = 0; i < N; i++)
      // {
      //     cout<<"arr = "<<soln_arr[i]<<endl;
      //     // soln_arr++;
      // }
      // cout<<"====================================="<<endl;
      // // std::copy(rand_arr, rand_arr + N, arr);
      // // avx2::SIMDSort(N, soln_arr);
      // for (size_t i = 0; i < N; i++)
      // {
      //     cout<<"arr = "<<soln_arr[i]<<endl;
      //     // soln_arr++;
      // }
      
      return 0;
  }
