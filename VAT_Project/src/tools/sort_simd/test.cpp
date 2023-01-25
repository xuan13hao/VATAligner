#include "./avx256/simd_sort.h"
#include "common.h"
#include <random>
#define LO -100
#define HI 100
#define NNUM 10
#define MAX std::numeric_limits<uint64_t>::max()
//g++ test.cpp common.cpp avx256/* -march=native -o a
using namespace std;
  // template <typename T>
  // static void RandGenInt(T* &arr, size_t N, T lo, T hi) 
  // {
  //   aligned_init<T>(arr, N);
  //   std::random_device rd;  //Will be used to obtain a seed for the random number engine
  //   std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  //   std::uniform_int_distribution<T> dis(lo, hi);
  //   for(size_t i = 0; i < N; i++) 
  //   {
  //     arr[i] = dis(gen);
  //   //   cout<<"arr = "<<arr[i]<<endl;
  //   }
  // }
  template <typename T>
  static void RandGenIntRecords(std::pair<T, T>* &arr, size_t N, T lo, T hi, unsigned int offset_start=0) {
    aligned_init<std::pair<T, T>>(arr, N);
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<T> dis(lo, hi);
    for(size_t i = 0; i < N; i++) {
      arr[i].first = dis(gen);
      arr[i].second = offset_start++;
    }
  }

  template <typename T>
  static void RandGenIntRecords(Tuple* &arr, size_t N, T lo, T hi, unsigned int offset_start=0) {
    aligned_init<Tuple>(arr, N);
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<T> dis(lo, hi);
    for(size_t i = 0; i < N; i++) {
      arr[i].key= dis(gen);
      arr[i].value = offset_start++;
    }
  }


int main(int argc, const char** argv) 
{
    using T = int64_t;
    size_t N = NNUM;
    int BLOCK_SIZE = 8;
    T lo = LO;
    T hi = HI;
    std::pair<T, T> *rand_arr;
    std::pair<T, T> *soln_arr;

    RandGenIntRecords(rand_arr, N, lo, hi);
    int p = N % 8;
    int n = N;
    if(p != 0)
    {   
       n = N + 8 - p  ;
    }
    vector<Tuple> tpl;
    for (size_t i = 0; i < n; i++)
    {
          if (i < N )
          {
            tpl.push_back(rand_arr[i]);
          }else
          {
            tpl.push_back(Tuple(MAX,MAX));
          }
          
    }
    for (size_t i = 0; i < tpl.size(); i++)
    {
        cout<<"arr "<<i<<"= "<<tpl[i].key<<":"<<tpl[i].value<<endl;
    }
    cout<<"====================================="<<endl;
    soln_arr = &tpl[0];
    aligned_init<Tuple>(soln_arr, tpl.size());
    std::copy(&tpl[0], &tpl[0] + tpl.size(), soln_arr);
    for (size_t i = 0; i < n; i++)
    {
        cout<<"arr "<<i<<"= "<<soln_arr[i].key<<":"<<soln_arr[i].value<<endl;
        // soln_arr++;
    }
    // cout<<"====================================="<<endl;
    // std::copy(rand_arr, rand_arr + N, arr);
    avx2::SIMDSort(n, soln_arr);
    for (size_t i = 0; i < tpl.size(); i++)
    {
        if ((soln_arr[i].key != MAX) && (soln_arr[i].value != MAX))
        {
          cout<<"arr "<<i<<"= "<<soln_arr[i].key<<":"<<soln_arr[i].value<<endl;
        }
    }
    delete rand_arr;
    delete soln_arr;
    return 0;
}