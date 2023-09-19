

#ifndef SIMD_SORT_H_GUARD
#define SIMD_SORT_H_GUARD

#include <vector>
#include <immintrin.h>
#include <assert.h>
#include <algorithm>

#ifndef AVX2_CAS_HPP
#define AVX2_CAS_HPP

#include <cstddef>
#include <immintrin.h>

#define AVX2_CAS_NO_BRANCHING 1
#define AVX2_CAS_ENABLE_SIMD 1

namespace avx2_cas // Compare and swap (cas).
{
    template <typename T, size_t const a0, size_t const b0> void cas (T *data);

    // Declared as a class to allow partial-specialization.
    template <typename T> class cas2 {public: template <size_t const a0, size_t const b0, size_t const a1, size_t const b1> static void cas (T *data);};
    template <typename T> class cas3 {public: template <size_t const a0, size_t const b0, size_t const a1, size_t const b1, size_t const a2, size_t const b2> static void cas (T *data);};
    template <typename T> class cas4 {public: template <size_t const a0, size_t const b0, size_t const a1, size_t const b1, size_t const a2, size_t const b2, size_t const a3, size_t const b3> static void cas (T *data);};

} // end namespace avx2_cas


template <typename T, size_t const a0, size_t const b0>
void
avx2_cas::cas(T * const data)
{
    static_assert(a0 < b0, "Function vector_sort::internal::cas requires i to be less than j.");
    T &d0 = data[a0];
    T &d1 = data[b0];
#if AVX2_CAS_NO_BRANCHING
    int const pi = ((int)(d0 < d1));
    int const qi = ((int)1) - pi;
    T const p = static_cast<T>(pi);
    T const q = static_cast<T>(qi);
    T const a = (p * d0) + (q * d1);
    T const b = (q * d0) + (p * d1);
    d0 = a;
    d1 = b;
#else
    if (d1 < d0)
    {
        std::swap(d0, d1);
    }
#endif
}

template<>
template <size_t const a0, size_t const b0, size_t const a1, size_t const b1>
void
avx2_cas::cas2<double>::cas(double * const data)
{
#if AVX2_CAS_ENABLE_SIMD
    static_assert(a0 < b0, "Function avx2_cas::cas2<double>::cas requires a0 to be less than b0.");
    static_assert(a1 < b1, "Function avx2_cas::cas2<double>::cas requires a1 to be less than b1.");

    // Load the 4 values into two vectors.
    __m128d A;
    __m128d B;
    double * const ap = reinterpret_cast<double *>(& A);
    double * const bp = reinterpret_cast<double *>(& B);
    ap[0] = data[a0];
    ap[1] = data[a1];
    bp[0] = data[b0];
    bp[1] = data[b1];

    // Compare the vectors.
    __m128d AA = _mm_min_pd(A, B);
    __m128d BB = _mm_max_pd(A, B);

    // Store the values back into the source array.
    double * const aap = reinterpret_cast<double *>(& AA);
    double * const bbp = reinterpret_cast<double *>(& BB);
    data[a0] = aap[0];
    data[a1] = aap[1];
    data[b0] = bbp[0];
    data[b1] = bbp[1];
#else
    avx2_cas::cas<double, a0, b0>(data);
    avx2_cas::cas<double, a1, b1>(data);
#endif
}

template<>
template <size_t const a0, size_t const b0, size_t const a1, size_t const b1>
void
avx2_cas::cas2<float>::cas(float * const data)
{
#if AVX2_CAS_ENABLE_SIMD
    static_assert(a0 < b0, "Function avx2_cas::cas2<float>::cas requires a0 to be less than b0.");
    static_assert(a1 < b1, "Function avx2_cas::cas2<float>::cas requires a1 to be less than b1.");

    // Load the 4 values into two vectors.
    __m128 A;
    __m128 B;
    float * const ap = reinterpret_cast<float *>(& A);
    float * const bp = reinterpret_cast<float *>(& B);
    ap[0] = data[a0];
    ap[1] = data[a1];
    ap[2] = 0.0F;
    ap[3] = 0.0F;
    bp[0] = data[b0];
    bp[1] = data[b1];
    bp[2] = 0.0F;
    bp[3] = 0.0F;

    // Compare the vectors.
    __m128 AA = _mm_min_ps(A, B);
    __m128 BB = _mm_max_ps(A, B);

    // Store the values back into the source array.
    float * const aap = reinterpret_cast<float *>(& AA);
    float * const bbp = reinterpret_cast<float *>(& BB);
    data[a0] = aap[0];
    data[a1] = aap[1];
    data[b0] = bbp[0];
    data[b1] = bbp[1];
#else
    avx2_cas::cas<float, a0, b0>(data);
    avx2_cas::cas<float, a1, b1>(data);
#endif
}

template<>
template <size_t const a0, size_t const b0, size_t const a1, size_t const b1>
void
avx2_cas::cas2<unsigned>::cas(unsigned * const data)
{
#if AVX2_CAS_ENABLE_SIMD
    static_assert(a0 < b0, "Function avx2_cas::cas2<unsigned>::cas requires a0 to be less than b0.");
    static_assert(a1 < b1, "Function avx2_cas::cas2<unsigned>::cas requires a1 to be less than b1.");

    // Load the 4 values into two vectors.
    __m128i A;
    __m128i B;
    unsigned * const ap = reinterpret_cast<unsigned *>(& A);
    unsigned * const bp = reinterpret_cast<unsigned *>(& B);
    ap[0] = data[a0];
    ap[1] = data[a1];
    ap[2] = 0;
    ap[3] = 0;
    bp[0] = data[b0];
    bp[1] = data[b1];
    bp[2] = 0;
    bp[3] = 0;

    // Compare the vectors.
    __m128i AA = _mm_min_epi32(A, B);
    __m128i BB = _mm_max_epi32(A, B);

    // Store the values back into the source array.
    unsigned * const aap = reinterpret_cast<unsigned *>(& AA);
    unsigned * const bbp = reinterpret_cast<unsigned *>(& BB);
    data[a0] = aap[0];
    data[a1] = aap[1];
    data[b0] = bbp[0];
    data[b1] = bbp[1];
#else
    avx2_cas::cas<unsigned, a0, b0>(data);
    avx2_cas::cas<unsigned, a1, b1>(data);
#endif
}

template<>
template <size_t const a0, size_t const b0, size_t const a1, size_t const b1, size_t const a2, size_t const b2>
void
avx2_cas::cas3<double>::cas(double * const data)
{
#if AVX2_CAS_ENABLE_SIMD
    static_assert(a0 < b0, "Function avx2_cas::cas2<double>::cas requires a0 to be less than b0.");
    static_assert(a1 < b1, "Function avx2_cas::cas2<double>::cas requires a1 to be less than b1.");

    // Load the 4 values into two vectors.
    __m256d A;
    __m256d B;
    double * const ap = reinterpret_cast<double *>(& A);
    double * const bp = reinterpret_cast<double *>(& B);
    ap[0] = data[a0];
    ap[1] = data[a1];
    ap[2] = data[a2];
    ap[3] = 0.0;
    bp[0] = data[b0];
    bp[1] = data[b1];
    bp[2] = data[b2];
    bp[3] = 0.0;

    // Compare the vectors.
    __m256d AA = _mm256_min_pd(A, B);
    __m256d BB = _mm256_max_pd(A, B);

    // Store the values back into the source array.
    double * const aap = reinterpret_cast<double *>(& AA);
    double * const bbp = reinterpret_cast<double *>(& BB);
    data[a0] = aap[0];
    data[a1] = aap[1];
    data[a2] = aap[2];
    data[b0] = bbp[0];
    data[b1] = bbp[1];
    data[b2] = bbp[2];
#else
    avx2_cas::cas<double, a0, b0>(data);
    avx2_cas::cas<double, a1, b1>(data);
    avx2_cas::cas<double, a2, b2>(data);
#endif
}

template<>
template <size_t const a0, size_t const b0, size_t const a1, size_t const b1, size_t const a2, size_t const b2>
void
avx2_cas::cas3<float>::cas(float * const data)
{
#if AVX2_CAS_ENABLE_SIMD
    static_assert(a0 < b0, "Function avx2_cas::cas3<float>::cas requires a0 to be less than b0.");
    static_assert(a1 < b1, "Function avx2_cas::cas3<float>::cas requires a1 to be less than b1.");
    static_assert(a2 < b2, "Function avx2_cas::cas3<float>::cas requires a2 to be less than b2.");

    // Load the 4 values into two vectors.
    __m128 A;
    __m128 B;
    float * const ap = reinterpret_cast<float *>(& A);
    float * const bp = reinterpret_cast<float *>(& B);
    ap[0] = data[a0];
    ap[1] = data[a1];
    ap[2] = data[a2];
    ap[3] = 0.0F;
    bp[0] = data[b0];
    bp[1] = data[b1];
    bp[2] = data[b2];
    bp[3] = 0.0F;

    // Compare the vectors.
    __m128 AA = _mm_min_ps(A, B);
    __m128 BB = _mm_max_ps(A, B);

    // Store the values back into the source array.
    float * const aap = reinterpret_cast<float *>(& AA);
    float * const bbp = reinterpret_cast<float *>(& BB);
    data[a0] = aap[0];
    data[a1] = aap[1];
    data[a2] = aap[2];
    data[b0] = bbp[0];
    data[b1] = bbp[1];
    data[b2] = bbp[2];
#else
    avx2_cas::cas<float, a0, b0>(data);
    avx2_cas::cas<float, a1, b1>(data);
    avx2_cas::cas<float, a2, b2>(data);
#endif
}

template<>
template <size_t const a0, size_t const b0, size_t const a1, size_t const b1, size_t const a2, size_t const b2>
void
avx2_cas::cas3<unsigned>::cas(unsigned * const data)
{
#if AVX2_CAS_ENABLE_SIMD
    static_assert(a0 < b0, "Function avx2_cas::cas2<unsigned>::cas requires a0 to be less than b0.");
    static_assert(a1 < b1, "Function avx2_cas::cas2<unsigned>::cas requires a1 to be less than b1.");
    static_assert(a2 < b2, "Function avx2_cas::cas3<unsigned>::cas requires a2 to be less than b2.");

    // Load the 4 values into two vectors.
    __m128i A;
    __m128i B;
    unsigned * const ap = reinterpret_cast<unsigned *>(& A);
    unsigned * const bp = reinterpret_cast<unsigned *>(& B);
    ap[0] = data[a0];
    ap[1] = data[a1];
    ap[2] = data[a2];
    ap[3] = 0;
    bp[0] = data[b0];
    bp[1] = data[b1];
    bp[2] = data[b2];
    bp[3] = 0;

    // Compare the vectors.
    __m128i AA = _mm_min_epi32(A, B);
    __m128i BB = _mm_max_epi32(A, B);

    // Store the values back into the source array.
    unsigned * const aap = reinterpret_cast<unsigned *>(& AA);
    unsigned * const bbp = reinterpret_cast<unsigned *>(& BB);
    data[a0] = aap[0];
    data[a1] = aap[1];
    data[a2] = aap[2];
    data[b0] = bbp[0];
    data[b1] = bbp[1];
    data[b2] = bbp[2];
#else
    avx2_cas::cas<unsigned, a0, b0>(data);
    avx2_cas::cas<unsigned, a1, b1>(data);
    avx2_cas::cas<unsigned, a2, b2>(data);
#endif
}

template<>
template <size_t const a0, size_t const b0, size_t const a1, size_t const b1, size_t const a2, size_t const b2, size_t const a3, size_t const b3>
void
avx2_cas::cas4<double>::cas(double * const data)
{
#if AVX2_CAS_ENABLE_SIMD
    static_assert(a0 < b0, "Function avx2_cas::cas4<double>::cas requires a0 to be less than b0.");
    static_assert(a1 < b1, "Function avx2_cas::cas4<double>::cas requires a1 to be less than b1.");
    static_assert(a2 < b2, "Function avx2_cas::cas4<double>::cas requires a2 to be less than b2.");
    static_assert(a3 < b3, "Function avx2_cas::cas4<double>::cas requires a3 to be less than b3.");

    // Load the 4 values into two vectors.
    __m256d A;
    __m256d B;
    double * const ap = reinterpret_cast<double *>(& A);
    double * const bp = reinterpret_cast<double *>(& B);
    ap[0] = data[a0];
    ap[1] = data[a1];
    ap[2] = data[a2];
    ap[3] = data[a3];
    bp[0] = data[b0];
    bp[1] = data[b1];
    bp[2] = data[b2];
    bp[3] = data[b3];

    // Compare the vectors.
    __m256d AA = _mm256_min_pd(A, B);
    __m256d BB = _mm256_max_pd(A, B);

    // Store the values back into the source array.
    double * const aap = reinterpret_cast<double *>(& AA);
    double * const bbp = reinterpret_cast<double *>(& BB);
    data[a0] = aap[0];
    data[a1] = aap[1];
    data[a2] = aap[2];
    data[a3] = aap[3];
    data[b0] = bbp[0];
    data[b1] = bbp[1];
    data[b2] = bbp[2];
    data[b3] = bbp[3];
#else
    avx2_cas::cas<double, a0, b0>(data);
    avx2_cas::cas<double, a1, b1>(data);
    avx2_cas::cas<double, a2, b2>(data);
    avx2_cas::cas<double, a3, b3>(data);
#endif
}

template<>
template <size_t const a0, size_t const b0, size_t const a1, size_t const b1, size_t const a2, size_t const b2, size_t const a3, size_t const b3>
void
avx2_cas::cas4<float>::cas(float * const data)
{
#if AVX2_CAS_ENABLE_SIMD
    static_assert(a0 < b0, "Function avx2_cas::cas4<double>::cas requires a0 to be less than b0.");
    static_assert(a1 < b1, "Function avx2_cas::cas4<double>::cas requires a1 to be less than b1.");
    static_assert(a2 < b2, "Function avx2_cas::cas4<double>::cas requires a2 to be less than b2.");
    static_assert(a3 < b3, "Function avx2_cas::cas4<double>::cas requires a3 to be less than b3.");

    // Load the 4 values into two vectors.
    __m128 A;
    __m128 B;
    float * const ap = reinterpret_cast<float *>(& A);
    float * const bp = reinterpret_cast<float *>(& B);
    ap[0] = data[a0];
    ap[1] = data[a1];
    ap[2] = data[a2];
    ap[3] = data[a3];
    bp[0] = data[b0];
    bp[1] = data[b1];
    bp[2] = data[b2];
    bp[3] = data[b3];

    // Compare the vectors.
    __m128 AA = _mm_min_ps(A, B);
    __m128 BB = _mm_max_ps(A, B);

    // Store the values back into the source array.
    float * const aap = reinterpret_cast<float *>(& AA);
    float * const bbp = reinterpret_cast<float *>(& BB);
    data[a0] = aap[0];
    data[a1] = aap[1];
    data[a2] = aap[2];
    data[a3] = aap[3];
    data[b0] = bbp[0];
    data[b1] = bbp[1];
    data[b2] = bbp[2];
    data[b3] = bbp[3];
#else
    avx2_cas::cas<float, a0, b0>(data);
    avx2_cas::cas<float, a1, b1>(data);
    avx2_cas::cas<float, a2, b2>(data);
    avx2_cas::cas<float, a3, b3>(data);
#endif
}

template<>
template <size_t const a0, size_t const b0, size_t const a1, size_t const b1, size_t const a2, size_t const b2, size_t const a3, size_t const b3>
void
avx2_cas::cas4<unsigned>::cas(unsigned * const data)
{
#if AVX2_CAS_ENABLE_SIMD
    static_assert(a0 < b0, "Function avx2_cas::cas4<unsigned>::cas requires a0 to be less than b0.");
    static_assert(a1 < b1, "Function avx2_cas::cas4<unsigned>::cas requires a1 to be less than b1.");
    static_assert(a2 < b2, "Function avx2_cas::cas4<unsigned>::cas requires a2 to be less than b2.");
    static_assert(a3 < b3, "Function avx2_cas::cas4<unsigned>::cas requires a3 to be less than b3.");

    // Load the 4 values into two vectors.
    __m128i A;
    __m128i B;
    unsigned * const ap = reinterpret_cast<unsigned *>(& A);
    unsigned * const bp = reinterpret_cast<unsigned *>(& B);
    ap[0] = data[a0];
    ap[1] = data[a1];
    ap[2] = data[a2];
    ap[3] = data[a3];
    bp[0] = data[b0];
    bp[1] = data[b1];
    bp[2] = data[b2];
    bp[3] = data[b3];

    // Compare the vectors.
    __m128i AA = _mm_min_epi32(A, B);
    __m128i BB = _mm_max_epi32(A, B);

    // Store the values back into the source array.
    unsigned * const aap = reinterpret_cast<unsigned *>(& AA);
    unsigned * const bbp = reinterpret_cast<unsigned *>(& BB);
    data[a0] = aap[0];
    data[a1] = aap[1];
    data[a2] = aap[2];
    data[a3] = aap[3];
    data[b0] = bbp[0];
    data[b1] = bbp[1];
    data[b2] = bbp[2];
    data[b3] = bbp[3];
#else
    avx2_cas::cas<unsigned, a0, b0>(data);
    avx2_cas::cas<unsigned, a1, b1>(data);
    avx2_cas::cas<unsigned, a2, b2>(data);
    avx2_cas::cas<unsigned, a3, b3>(data);
#endif
}

#endif // AVX2_CAS_HPP

 namespace vector_sort{ namespace internal{ unsigned long switchPoint = 32; } } 

/* 
        |           |   _)                               |   
   __|  __|   _` |  __|  |   __|       __|   _ \    __|  __|
 \__ \  |    (   |  |    |  (        \__ \  (   |  |     |   
 ____/ \__| \__,_| \__| _| \___|     ____/ \___/  _|    \__|   
                   Copyright (c) 2018 M. Welsch <michael@welsch.one>
       Copyright (c) 2018 S. D. Adams <s.d.adams.software@gmail.com> 
*/

#ifndef SORTING_NETWORK_HPP
#define SORTING_NETWORK_HPP

// This file is automatically generated by dev/avx_static_sort_generator.cpp.  Do not edit manually, changes will be lost.

#ifndef SIMD_SORT_H_GUARD
#include "avx2_cas.cpp"

#endif
namespace vector_sort
{
    template <typename T> void sort(T *data, size_t n);

    namespace internal
    {
        template <size_t const n> class net {public: template <typename T> static void sort(T *);};

        template <typename T> 
        inline void static_sort(T * const data, size_t const n); 
    }
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<2>::sort(T * const data)
{
    avx2_cas::cas<T, 0, 1>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<3>::sort(T * const data)
{
    avx2_cas::cas<T, 1, 2>(data);
    avx2_cas::cas<T, 0, 2>(data);
    avx2_cas::cas<T, 0, 1>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<4>::sort(T * const data)
{
    avx2_cas::cas2<T>::template cas<0, 1, 2, 3>(data);
    avx2_cas::cas2<T>::template cas<0, 2, 1, 3>(data);
    avx2_cas::cas<T, 1, 2>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<5>::sort(T * const data)
{
    avx2_cas::cas2<T>::template cas<0, 1, 3, 4>(data);
    avx2_cas::cas<T, 2, 4>(data);
    avx2_cas::cas<T, 2, 3>(data);
    avx2_cas::cas<T, 0, 3>(data);
    avx2_cas::cas2<T>::template cas<0, 2, 1, 4>(data);
    avx2_cas::cas<T, 1, 3>(data);
    avx2_cas::cas<T, 1, 2>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<6>::sort(T * const data)
{
    avx2_cas::cas<T, 1, 2>(data);
    avx2_cas::cas<T, 0, 2>(data);
    avx2_cas::cas2<T>::template cas<0, 1, 4, 5>(data);
    avx2_cas::cas<T, 3, 5>(data);
    avx2_cas::cas<T, 3, 4>(data);
    avx2_cas::cas3<T>::template cas<0, 3, 1, 4, 2, 5>(data);
    avx2_cas::cas2<T>::template cas<2, 4, 1, 3>(data);
    avx2_cas::cas<T, 2, 3>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<7>::sort(T * const data)
{
    avx2_cas::cas<T, 1, 2>(data);
    avx2_cas::cas<T, 0, 2>(data);
    avx2_cas::cas3<T>::template cas<0, 1, 3, 4, 5, 6>(data);
    avx2_cas::cas2<T>::template cas<3, 5, 4, 6>(data);
    avx2_cas::cas<T, 4, 5>(data);
    avx2_cas::cas<T, 0, 4>(data);
    avx2_cas::cas3<T>::template cas<0, 3, 1, 5, 2, 6>(data);
    avx2_cas::cas2<T>::template cas<2, 5, 1, 3>(data);
    avx2_cas::cas<T, 2, 4>(data);
    avx2_cas::cas<T, 2, 3>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<8>::sort(T * const data)
{
    avx2_cas::cas2<T>::template cas<0, 1, 2, 3>(data);
    avx2_cas::cas2<T>::template cas<0, 2, 1, 3>(data);
    avx2_cas::cas3<T>::template cas<1, 2, 4, 5, 6, 7>(data);
    avx2_cas::cas2<T>::template cas<4, 6, 5, 7>(data);
    avx2_cas::cas2<T>::template cas<5, 6, 0, 4>(data);
    avx2_cas::cas<T, 1, 5>(data);
    avx2_cas::cas3<T>::template cas<1, 4, 2, 6, 3, 7>(data);
    avx2_cas::cas2<T>::template cas<3, 6, 2, 4>(data);
    avx2_cas::cas<T, 3, 5>(data);
    avx2_cas::cas<T, 3, 4>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<9>::sort(T * const data)
{
    avx2_cas::cas2<T>::template cas<0, 1, 2, 3>(data);
    avx2_cas::cas2<T>::template cas<0, 2, 1, 3>(data);
    avx2_cas::cas3<T>::template cas<1, 2, 4, 5, 7, 8>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas<T, 6, 7>(data);
    avx2_cas::cas<T, 4, 7>(data);
    avx2_cas::cas2<T>::template cas<4, 6, 5, 8>(data);
    avx2_cas::cas<T, 5, 7>(data);
    avx2_cas::cas<T, 5, 6>(data);
    avx2_cas::cas<T, 0, 5>(data);
    avx2_cas::cas2<T>::template cas<0, 4, 1, 6>(data);
    avx2_cas::cas<T, 1, 5>(data);
    avx2_cas::cas3<T>::template cas<1, 4, 2, 7, 3, 8>(data);
    avx2_cas::cas2<T>::template cas<3, 7, 2, 5>(data);
    avx2_cas::cas2<T>::template cas<2, 4, 3, 6>(data);
    avx2_cas::cas<T, 3, 5>(data);
    avx2_cas::cas<T, 3, 4>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<10>::sort(T * const data)
{
    avx2_cas::cas2<T>::template cas<0, 1, 3, 4>(data);
    avx2_cas::cas<T, 2, 4>(data);
    avx2_cas::cas<T, 2, 3>(data);
    avx2_cas::cas<T, 0, 3>(data);
    avx2_cas::cas2<T>::template cas<0, 2, 1, 4>(data);
    avx2_cas::cas<T, 1, 3>(data);
    avx2_cas::cas3<T>::template cas<1, 2, 5, 6, 8, 9>(data);
    avx2_cas::cas<T, 7, 9>(data);
    avx2_cas::cas<T, 7, 8>(data);
    avx2_cas::cas<T, 5, 8>(data);
    avx2_cas::cas2<T>::template cas<5, 7, 6, 9>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas2<T>::template cas<6, 7, 0, 5>(data);
    avx2_cas::cas<T, 1, 6>(data);
    avx2_cas::cas4<T>::template cas<1, 5, 2, 7, 3, 8, 4, 9>(data);
    avx2_cas::cas2<T>::template cas<4, 8, 3, 7>(data);
    avx2_cas::cas3<T>::template cas<4, 7, 2, 5, 3, 6>(data);
    avx2_cas::cas2<T>::template cas<4, 6, 3, 5>(data);
    avx2_cas::cas<T, 4, 5>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<11>::sort(T * const data)
{
    avx2_cas::cas2<T>::template cas<0, 1, 3, 4>(data);
    avx2_cas::cas<T, 2, 4>(data);
    avx2_cas::cas<T, 2, 3>(data);
    avx2_cas::cas<T, 0, 3>(data);
    avx2_cas::cas2<T>::template cas<0, 2, 1, 4>(data);
    avx2_cas::cas<T, 1, 3>(data);
    avx2_cas::cas2<T>::template cas<1, 2, 6, 7>(data);
    avx2_cas::cas<T, 5, 7>(data);
    avx2_cas::cas2<T>::template cas<5, 6, 9, 10>(data);
    avx2_cas::cas<T, 8, 10>(data);
    avx2_cas::cas<T, 8, 9>(data);
    avx2_cas::cas3<T>::template cas<5, 8, 6, 9, 7, 10>(data);
    avx2_cas::cas2<T>::template cas<7, 9, 6, 8>(data);
    avx2_cas::cas2<T>::template cas<7, 8, 0, 6>(data);
    avx2_cas::cas2<T>::template cas<0, 5, 1, 7>(data);
    avx2_cas::cas<T, 1, 6>(data);
    avx2_cas::cas4<T>::template cas<1, 5, 2, 8, 3, 9, 4, 10>(data);
    avx2_cas::cas2<T>::template cas<4, 9, 3, 8>(data);
    avx2_cas::cas3<T>::template cas<4, 8, 2, 5, 3, 6>(data);
    avx2_cas::cas<T, 4, 7>(data);
    avx2_cas::cas2<T>::template cas<4, 6, 3, 5>(data);
    avx2_cas::cas<T, 4, 5>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<12>::sort(T * const data)
{
    avx2_cas::cas<T, 1, 2>(data);
    avx2_cas::cas<T, 0, 2>(data);
    avx2_cas::cas2<T>::template cas<0, 1, 4, 5>(data);
    avx2_cas::cas<T, 3, 5>(data);
    avx2_cas::cas<T, 3, 4>(data);
    avx2_cas::cas3<T>::template cas<0, 3, 1, 4, 2, 5>(data);
    avx2_cas::cas2<T>::template cas<2, 4, 1, 3>(data);
    avx2_cas::cas2<T>::template cas<2, 3, 7, 8>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas2<T>::template cas<6, 7, 10, 11>(data);
    avx2_cas::cas<T, 9, 11>(data);
    avx2_cas::cas<T, 9, 10>(data);
    avx2_cas::cas3<T>::template cas<6, 9, 7, 10, 8, 11>(data);
    avx2_cas::cas2<T>::template cas<8, 10, 7, 9>(data);
    avx2_cas::cas3<T>::template cas<8, 9, 0, 6, 1, 7>(data);
    avx2_cas::cas<T, 2, 8>(data);
    avx2_cas::cas2<T>::template cas<2, 7, 1, 6>(data);
    avx2_cas::cas4<T>::template cas<2, 6, 3, 9, 4, 10, 5, 11>(data);
    avx2_cas::cas2<T>::template cas<5, 10, 4, 9>(data);
    avx2_cas::cas3<T>::template cas<5, 9, 3, 6, 4, 7>(data);
    avx2_cas::cas<T, 5, 8>(data);
    avx2_cas::cas2<T>::template cas<5, 7, 4, 6>(data);
    avx2_cas::cas<T, 5, 6>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<13>::sort(T * const data)
{
    avx2_cas::cas<T, 1, 2>(data);
    avx2_cas::cas<T, 0, 2>(data);
    avx2_cas::cas2<T>::template cas<0, 1, 4, 5>(data);
    avx2_cas::cas<T, 3, 5>(data);
    avx2_cas::cas<T, 3, 4>(data);
    avx2_cas::cas3<T>::template cas<0, 3, 1, 4, 2, 5>(data);
    avx2_cas::cas2<T>::template cas<2, 4, 1, 3>(data);
    avx2_cas::cas2<T>::template cas<2, 3, 7, 8>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas3<T>::template cas<6, 7, 9, 10, 11, 12>(data);
    avx2_cas::cas2<T>::template cas<9, 11, 10, 12>(data);
    avx2_cas::cas<T, 10, 11>(data);
    avx2_cas::cas<T, 6, 10>(data);
    avx2_cas::cas3<T>::template cas<6, 9, 7, 11, 8, 12>(data);
    avx2_cas::cas2<T>::template cas<8, 11, 7, 9>(data);
    avx2_cas::cas<T, 8, 10>(data);
    avx2_cas::cas2<T>::template cas<8, 9, 0, 7>(data);
    avx2_cas::cas3<T>::template cas<0, 6, 1, 8, 2, 9>(data);
    avx2_cas::cas2<T>::template cas<2, 8, 1, 6>(data);
    avx2_cas::cas<T, 2, 7>(data);
    avx2_cas::cas4<T>::template cas<2, 6, 3, 10, 4, 11, 5, 12>(data);
    avx2_cas::cas2<T>::template cas<5, 11, 4, 10>(data);
    avx2_cas::cas2<T>::template cas<5, 10, 3, 7>(data);
    avx2_cas::cas3<T>::template cas<3, 6, 4, 8, 5, 9>(data);
    avx2_cas::cas2<T>::template cas<5, 8, 4, 6>(data);
    avx2_cas::cas<T, 5, 7>(data);
    avx2_cas::cas<T, 5, 6>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<14>::sort(T * const data)
{
    avx2_cas::cas<T, 1, 2>(data);
    avx2_cas::cas<T, 0, 2>(data);
    avx2_cas::cas3<T>::template cas<0, 1, 3, 4, 5, 6>(data);
    avx2_cas::cas2<T>::template cas<3, 5, 4, 6>(data);
    avx2_cas::cas<T, 4, 5>(data);
    avx2_cas::cas<T, 0, 4>(data);
    avx2_cas::cas3<T>::template cas<0, 3, 1, 5, 2, 6>(data);
    avx2_cas::cas2<T>::template cas<2, 5, 1, 3>(data);
    avx2_cas::cas<T, 2, 4>(data);
    avx2_cas::cas2<T>::template cas<2, 3, 8, 9>(data);
    avx2_cas::cas<T, 7, 9>(data);
    avx2_cas::cas3<T>::template cas<7, 8, 10, 11, 12, 13>(data);
    avx2_cas::cas2<T>::template cas<10, 12, 11, 13>(data);
    avx2_cas::cas<T, 11, 12>(data);
    avx2_cas::cas<T, 7, 11>(data);
    avx2_cas::cas3<T>::template cas<7, 10, 8, 12, 9, 13>(data);
    avx2_cas::cas2<T>::template cas<9, 12, 8, 10>(data);
    avx2_cas::cas<T, 9, 11>(data);
    avx2_cas::cas3<T>::template cas<9, 10, 0, 7, 1, 8>(data);
    avx2_cas::cas<T, 2, 9>(data);
    avx2_cas::cas2<T>::template cas<2, 8, 1, 7>(data);
    avx2_cas::cas3<T>::template cas<2, 7, 3, 10, 4, 11>(data);
    avx2_cas::cas3<T>::template cas<4, 10, 5, 12, 6, 13>(data);
    avx2_cas::cas2<T>::template cas<6, 12, 5, 10>(data);
    avx2_cas::cas<T, 6, 11>(data);
    avx2_cas::cas3<T>::template cas<6, 10, 3, 7, 4, 8>(data);
    avx2_cas::cas2<T>::template cas<4, 7, 5, 9>(data);
    avx2_cas::cas2<T>::template cas<6, 9, 5, 7>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas<T, 6, 7>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<15>::sort(T * const data)
{
    avx2_cas::cas<T, 1, 2>(data);
    avx2_cas::cas<T, 0, 2>(data);
    avx2_cas::cas3<T>::template cas<0, 1, 3, 4, 5, 6>(data);
    avx2_cas::cas2<T>::template cas<3, 5, 4, 6>(data);
    avx2_cas::cas<T, 4, 5>(data);
    avx2_cas::cas<T, 0, 4>(data);
    avx2_cas::cas3<T>::template cas<0, 3, 1, 5, 2, 6>(data);
    avx2_cas::cas2<T>::template cas<2, 5, 1, 3>(data);
    avx2_cas::cas<T, 2, 4>(data);
    avx2_cas::cas3<T>::template cas<2, 3, 7, 8, 9, 10>(data);
    avx2_cas::cas2<T>::template cas<7, 9, 8, 10>(data);
    avx2_cas::cas3<T>::template cas<8, 9, 11, 12, 13, 14>(data);
    avx2_cas::cas2<T>::template cas<11, 13, 12, 14>(data);
    avx2_cas::cas2<T>::template cas<12, 13, 7, 11>(data);
    avx2_cas::cas<T, 8, 12>(data);
    avx2_cas::cas3<T>::template cas<8, 11, 9, 13, 10, 14>(data);
    avx2_cas::cas2<T>::template cas<10, 13, 9, 11>(data);
    avx2_cas::cas<T, 10, 12>(data);
    avx2_cas::cas2<T>::template cas<10, 11, 0, 8>(data);
    avx2_cas::cas3<T>::template cas<0, 7, 1, 9, 2, 10>(data);
    avx2_cas::cas2<T>::template cas<2, 9, 1, 7>(data);
    avx2_cas::cas<T, 2, 8>(data);
    avx2_cas::cas3<T>::template cas<2, 7, 3, 11, 4, 12>(data);
    avx2_cas::cas3<T>::template cas<4, 11, 5, 13, 6, 14>(data);
    avx2_cas::cas2<T>::template cas<6, 13, 5, 11>(data);
    avx2_cas::cas<T, 6, 12>(data);
    avx2_cas::cas3<T>::template cas<6, 11, 3, 7, 4, 8>(data);
    avx2_cas::cas3<T>::template cas<4, 7, 5, 9, 6, 10>(data);
    avx2_cas::cas2<T>::template cas<6, 9, 5, 7>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas<T, 6, 7>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<16>::sort(T * const data)
{
    avx2_cas::cas2<T>::template cas<0, 1, 2, 3>(data);
    avx2_cas::cas2<T>::template cas<0, 2, 1, 3>(data);
    avx2_cas::cas3<T>::template cas<1, 2, 4, 5, 6, 7>(data);
    avx2_cas::cas2<T>::template cas<4, 6, 5, 7>(data);
    avx2_cas::cas2<T>::template cas<5, 6, 0, 4>(data);
    avx2_cas::cas<T, 1, 5>(data);
    avx2_cas::cas3<T>::template cas<1, 4, 2, 6, 3, 7>(data);
    avx2_cas::cas2<T>::template cas<3, 6, 2, 4>(data);
    avx2_cas::cas<T, 3, 5>(data);
    avx2_cas::cas3<T>::template cas<3, 4, 8, 9, 10, 11>(data);
    avx2_cas::cas2<T>::template cas<8, 10, 9, 11>(data);
    avx2_cas::cas3<T>::template cas<9, 10, 12, 13, 14, 15>(data);
    avx2_cas::cas2<T>::template cas<12, 14, 13, 15>(data);
    avx2_cas::cas2<T>::template cas<13, 14, 8, 12>(data);
    avx2_cas::cas<T, 9, 13>(data);
    avx2_cas::cas3<T>::template cas<9, 12, 10, 14, 11, 15>(data);
    avx2_cas::cas2<T>::template cas<11, 14, 10, 12>(data);
    avx2_cas::cas<T, 11, 13>(data);
    avx2_cas::cas3<T>::template cas<11, 12, 0, 8, 1, 9>(data);
    avx2_cas::cas3<T>::template cas<1, 8, 2, 10, 3, 11>(data);
    avx2_cas::cas2<T>::template cas<3, 10, 2, 8>(data);
    avx2_cas::cas<T, 3, 9>(data);
    avx2_cas::cas3<T>::template cas<3, 8, 4, 12, 5, 13>(data);
    avx2_cas::cas3<T>::template cas<5, 12, 6, 14, 7, 15>(data);
    avx2_cas::cas2<T>::template cas<7, 14, 6, 12>(data);
    avx2_cas::cas<T, 7, 13>(data);
    avx2_cas::cas3<T>::template cas<7, 12, 4, 8, 5, 9>(data);
    avx2_cas::cas3<T>::template cas<5, 8, 6, 10, 7, 11>(data);
    avx2_cas::cas2<T>::template cas<7, 10, 6, 8>(data);
    avx2_cas::cas<T, 7, 9>(data);
    avx2_cas::cas<T, 7, 8>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<17>::sort(T * const data)
{
    avx2_cas::cas2<T>::template cas<0, 1, 2, 3>(data);
    avx2_cas::cas2<T>::template cas<0, 2, 1, 3>(data);
    avx2_cas::cas3<T>::template cas<1, 2, 4, 5, 6, 7>(data);
    avx2_cas::cas2<T>::template cas<4, 6, 5, 7>(data);
    avx2_cas::cas2<T>::template cas<5, 6, 0, 4>(data);
    avx2_cas::cas<T, 1, 5>(data);
    avx2_cas::cas3<T>::template cas<1, 4, 2, 6, 3, 7>(data);
    avx2_cas::cas2<T>::template cas<3, 6, 2, 4>(data);
    avx2_cas::cas<T, 3, 5>(data);
    avx2_cas::cas3<T>::template cas<3, 4, 8, 9, 10, 11>(data);
    avx2_cas::cas2<T>::template cas<8, 10, 9, 11>(data);
    avx2_cas::cas3<T>::template cas<9, 10, 12, 13, 15, 16>(data);
    avx2_cas::cas<T, 14, 16>(data);
    avx2_cas::cas<T, 14, 15>(data);
    avx2_cas::cas<T, 12, 15>(data);
    avx2_cas::cas2<T>::template cas<12, 14, 13, 16>(data);
    avx2_cas::cas<T, 13, 15>(data);
    avx2_cas::cas<T, 13, 14>(data);
    avx2_cas::cas<T, 8, 13>(data);
    avx2_cas::cas2<T>::template cas<8, 12, 9, 14>(data);
    avx2_cas::cas<T, 9, 13>(data);
    avx2_cas::cas3<T>::template cas<9, 12, 10, 15, 11, 16>(data);
    avx2_cas::cas2<T>::template cas<11, 15, 10, 13>(data);
    avx2_cas::cas2<T>::template cas<10, 12, 11, 14>(data);
    avx2_cas::cas<T, 11, 13>(data);
    avx2_cas::cas2<T>::template cas<11, 12, 0, 9>(data);
    avx2_cas::cas2<T>::template cas<0, 8, 1, 10>(data);
    avx2_cas::cas<T, 1, 9>(data);
    avx2_cas::cas3<T>::template cas<1, 8, 2, 11, 3, 12>(data);
    avx2_cas::cas2<T>::template cas<3, 11, 2, 9>(data);
    avx2_cas::cas2<T>::template cas<2, 8, 3, 10>(data);
    avx2_cas::cas<T, 3, 9>(data);
    avx2_cas::cas3<T>::template cas<3, 8, 4, 13, 5, 14>(data);
    avx2_cas::cas3<T>::template cas<5, 13, 6, 15, 7, 16>(data);
    avx2_cas::cas2<T>::template cas<7, 15, 6, 13>(data);
    avx2_cas::cas<T, 7, 14>(data);
    avx2_cas::cas2<T>::template cas<7, 13, 4, 9>(data);
    avx2_cas::cas2<T>::template cas<4, 8, 5, 10>(data);
    avx2_cas::cas<T, 5, 9>(data);
    avx2_cas::cas3<T>::template cas<5, 8, 6, 11, 7, 12>(data);
    avx2_cas::cas2<T>::template cas<7, 11, 6, 9>(data);
    avx2_cas::cas2<T>::template cas<6, 8, 7, 10>(data);
    avx2_cas::cas<T, 7, 9>(data);
    avx2_cas::cas<T, 7, 8>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<18>::sort(T * const data)
{
    avx2_cas::cas2<T>::template cas<0, 1, 2, 3>(data);
    avx2_cas::cas2<T>::template cas<0, 2, 1, 3>(data);
    avx2_cas::cas3<T>::template cas<1, 2, 4, 5, 7, 8>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas<T, 6, 7>(data);
    avx2_cas::cas<T, 4, 7>(data);
    avx2_cas::cas2<T>::template cas<4, 6, 5, 8>(data);
    avx2_cas::cas<T, 5, 7>(data);
    avx2_cas::cas<T, 5, 6>(data);
    avx2_cas::cas<T, 0, 5>(data);
    avx2_cas::cas2<T>::template cas<0, 4, 1, 6>(data);
    avx2_cas::cas<T, 1, 5>(data);
    avx2_cas::cas3<T>::template cas<1, 4, 2, 7, 3, 8>(data);
    avx2_cas::cas2<T>::template cas<3, 7, 2, 5>(data);
    avx2_cas::cas2<T>::template cas<2, 4, 3, 6>(data);
    avx2_cas::cas<T, 3, 5>(data);
    avx2_cas::cas3<T>::template cas<3, 4, 9, 10, 11, 12>(data);
    avx2_cas::cas2<T>::template cas<9, 11, 10, 12>(data);
    avx2_cas::cas3<T>::template cas<10, 11, 13, 14, 16, 17>(data);
    avx2_cas::cas<T, 15, 17>(data);
    avx2_cas::cas<T, 15, 16>(data);
    avx2_cas::cas<T, 13, 16>(data);
    avx2_cas::cas2<T>::template cas<13, 15, 14, 17>(data);
    avx2_cas::cas<T, 14, 16>(data);
    avx2_cas::cas<T, 14, 15>(data);
    avx2_cas::cas<T, 9, 14>(data);
    avx2_cas::cas2<T>::template cas<9, 13, 10, 15>(data);
    avx2_cas::cas<T, 10, 14>(data);
    avx2_cas::cas3<T>::template cas<10, 13, 11, 16, 12, 17>(data);
    avx2_cas::cas2<T>::template cas<12, 16, 11, 14>(data);
    avx2_cas::cas2<T>::template cas<11, 13, 12, 15>(data);
    avx2_cas::cas<T, 12, 14>(data);
    avx2_cas::cas3<T>::template cas<12, 13, 0, 9, 1, 10>(data);
    avx2_cas::cas3<T>::template cas<1, 9, 2, 11, 3, 12>(data);
    avx2_cas::cas2<T>::template cas<3, 11, 2, 9>(data);
    avx2_cas::cas<T, 3, 10>(data);
    avx2_cas::cas3<T>::template cas<3, 9, 4, 13, 5, 14>(data);
    avx2_cas::cas4<T>::template cas<5, 13, 6, 15, 7, 16, 8, 17>(data);
    avx2_cas::cas2<T>::template cas<8, 16, 7, 15>(data);
    avx2_cas::cas3<T>::template cas<8, 15, 6, 13, 7, 14>(data);
    avx2_cas::cas2<T>::template cas<8, 14, 7, 13>(data);
    avx2_cas::cas3<T>::template cas<8, 13, 4, 9, 5, 10>(data);
    avx2_cas::cas3<T>::template cas<5, 9, 6, 11, 7, 12>(data);
    avx2_cas::cas2<T>::template cas<8, 12, 7, 11>(data);
    avx2_cas::cas3<T>::template cas<8, 11, 6, 9, 7, 10>(data);
    avx2_cas::cas2<T>::template cas<8, 10, 7, 9>(data);
    avx2_cas::cas<T, 8, 9>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<19>::sort(T * const data)
{
    avx2_cas::cas2<T>::template cas<0, 1, 2, 3>(data);
    avx2_cas::cas2<T>::template cas<0, 2, 1, 3>(data);
    avx2_cas::cas3<T>::template cas<1, 2, 4, 5, 7, 8>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas<T, 6, 7>(data);
    avx2_cas::cas<T, 4, 7>(data);
    avx2_cas::cas2<T>::template cas<4, 6, 5, 8>(data);
    avx2_cas::cas<T, 5, 7>(data);
    avx2_cas::cas<T, 5, 6>(data);
    avx2_cas::cas<T, 0, 5>(data);
    avx2_cas::cas2<T>::template cas<0, 4, 1, 6>(data);
    avx2_cas::cas<T, 1, 5>(data);
    avx2_cas::cas3<T>::template cas<1, 4, 2, 7, 3, 8>(data);
    avx2_cas::cas2<T>::template cas<3, 7, 2, 5>(data);
    avx2_cas::cas2<T>::template cas<2, 4, 3, 6>(data);
    avx2_cas::cas<T, 3, 5>(data);
    avx2_cas::cas3<T>::template cas<3, 4, 9, 10, 12, 13>(data);
    avx2_cas::cas<T, 11, 13>(data);
    avx2_cas::cas<T, 11, 12>(data);
    avx2_cas::cas<T, 9, 12>(data);
    avx2_cas::cas2<T>::template cas<9, 11, 10, 13>(data);
    avx2_cas::cas<T, 10, 12>(data);
    avx2_cas::cas3<T>::template cas<10, 11, 14, 15, 17, 18>(data);
    avx2_cas::cas<T, 16, 18>(data);
    avx2_cas::cas<T, 16, 17>(data);
    avx2_cas::cas<T, 14, 17>(data);
    avx2_cas::cas2<T>::template cas<14, 16, 15, 18>(data);
    avx2_cas::cas<T, 15, 17>(data);
    avx2_cas::cas2<T>::template cas<15, 16, 9, 14>(data);
    avx2_cas::cas<T, 10, 15>(data);
    avx2_cas::cas4<T>::template cas<10, 14, 11, 16, 12, 17, 13, 18>(data);
    avx2_cas::cas2<T>::template cas<13, 17, 12, 16>(data);
    avx2_cas::cas3<T>::template cas<13, 16, 11, 14, 12, 15>(data);
    avx2_cas::cas2<T>::template cas<13, 15, 12, 14>(data);
    avx2_cas::cas2<T>::template cas<13, 14, 0, 10>(data);
    avx2_cas::cas2<T>::template cas<0, 9, 1, 11>(data);
    avx2_cas::cas<T, 1, 10>(data);
    avx2_cas::cas3<T>::template cas<1, 9, 2, 12, 3, 13>(data);
    avx2_cas::cas2<T>::template cas<3, 12, 2, 10>(data);
    avx2_cas::cas2<T>::template cas<2, 9, 3, 11>(data);
    avx2_cas::cas<T, 3, 10>(data);
    avx2_cas::cas3<T>::template cas<3, 9, 4, 14, 5, 15>(data);
    avx2_cas::cas4<T>::template cas<5, 14, 6, 16, 7, 17, 8, 18>(data);
    avx2_cas::cas2<T>::template cas<8, 17, 7, 16>(data);
    avx2_cas::cas3<T>::template cas<8, 16, 6, 14, 7, 15>(data);
    avx2_cas::cas2<T>::template cas<8, 15, 7, 14>(data);
    avx2_cas::cas3<T>::template cas<8, 14, 4, 9, 5, 10>(data);
    avx2_cas::cas4<T>::template cas<5, 9, 6, 11, 7, 12, 8, 13>(data);
    avx2_cas::cas2<T>::template cas<8, 12, 7, 11>(data);
    avx2_cas::cas3<T>::template cas<8, 11, 6, 9, 7, 10>(data);
    avx2_cas::cas2<T>::template cas<8, 10, 7, 9>(data);
    avx2_cas::cas<T, 8, 9>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<20>::sort(T * const data)
{
    avx2_cas::cas2<T>::template cas<0, 1, 3, 4>(data);
    avx2_cas::cas<T, 2, 4>(data);
    avx2_cas::cas<T, 2, 3>(data);
    avx2_cas::cas<T, 0, 3>(data);
    avx2_cas::cas2<T>::template cas<0, 2, 1, 4>(data);
    avx2_cas::cas<T, 1, 3>(data);
    avx2_cas::cas3<T>::template cas<1, 2, 5, 6, 8, 9>(data);
    avx2_cas::cas<T, 7, 9>(data);
    avx2_cas::cas<T, 7, 8>(data);
    avx2_cas::cas<T, 5, 8>(data);
    avx2_cas::cas2<T>::template cas<5, 7, 6, 9>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas2<T>::template cas<6, 7, 0, 5>(data);
    avx2_cas::cas<T, 1, 6>(data);
    avx2_cas::cas4<T>::template cas<1, 5, 2, 7, 3, 8, 4, 9>(data);
    avx2_cas::cas2<T>::template cas<4, 8, 3, 7>(data);
    avx2_cas::cas3<T>::template cas<4, 7, 2, 5, 3, 6>(data);
    avx2_cas::cas2<T>::template cas<4, 6, 3, 5>(data);
    avx2_cas::cas3<T>::template cas<4, 5, 10, 11, 13, 14>(data);
    avx2_cas::cas<T, 12, 14>(data);
    avx2_cas::cas<T, 12, 13>(data);
    avx2_cas::cas<T, 10, 13>(data);
    avx2_cas::cas2<T>::template cas<10, 12, 11, 14>(data);
    avx2_cas::cas<T, 11, 13>(data);
    avx2_cas::cas3<T>::template cas<11, 12, 15, 16, 18, 19>(data);
    avx2_cas::cas<T, 17, 19>(data);
    avx2_cas::cas<T, 17, 18>(data);
    avx2_cas::cas<T, 15, 18>(data);
    avx2_cas::cas2<T>::template cas<15, 17, 16, 19>(data);
    avx2_cas::cas<T, 16, 18>(data);
    avx2_cas::cas2<T>::template cas<16, 17, 10, 15>(data);
    avx2_cas::cas<T, 11, 16>(data);
    avx2_cas::cas4<T>::template cas<11, 15, 12, 17, 13, 18, 14, 19>(data);
    avx2_cas::cas2<T>::template cas<14, 18, 13, 17>(data);
    avx2_cas::cas3<T>::template cas<14, 17, 12, 15, 13, 16>(data);
    avx2_cas::cas2<T>::template cas<14, 16, 13, 15>(data);
    avx2_cas::cas3<T>::template cas<14, 15, 0, 10, 1, 11>(data);
    avx2_cas::cas4<T>::template cas<1, 10, 2, 12, 3, 13, 4, 14>(data);
    avx2_cas::cas2<T>::template cas<4, 13, 3, 12>(data);
    avx2_cas::cas3<T>::template cas<4, 12, 2, 10, 3, 11>(data);
    avx2_cas::cas2<T>::template cas<4, 11, 3, 10>(data);
    avx2_cas::cas3<T>::template cas<4, 10, 5, 15, 6, 16>(data);
    avx2_cas::cas4<T>::template cas<6, 15, 7, 17, 8, 18, 9, 19>(data);
    avx2_cas::cas2<T>::template cas<9, 18, 8, 17>(data);
    avx2_cas::cas3<T>::template cas<9, 17, 7, 15, 8, 16>(data);
    avx2_cas::cas2<T>::template cas<9, 16, 8, 15>(data);
    avx2_cas::cas3<T>::template cas<9, 15, 5, 10, 6, 11>(data);
    avx2_cas::cas4<T>::template cas<6, 10, 7, 12, 8, 13, 9, 14>(data);
    avx2_cas::cas2<T>::template cas<9, 13, 8, 12>(data);
    avx2_cas::cas3<T>::template cas<9, 12, 7, 10, 8, 11>(data);
    avx2_cas::cas2<T>::template cas<9, 11, 8, 10>(data);
    avx2_cas::cas<T, 9, 10>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<21>::sort(T * const data)
{
    avx2_cas::cas2<T>::template cas<0, 1, 3, 4>(data);
    avx2_cas::cas<T, 2, 4>(data);
    avx2_cas::cas<T, 2, 3>(data);
    avx2_cas::cas<T, 0, 3>(data);
    avx2_cas::cas2<T>::template cas<0, 2, 1, 4>(data);
    avx2_cas::cas<T, 1, 3>(data);
    avx2_cas::cas3<T>::template cas<1, 2, 5, 6, 8, 9>(data);
    avx2_cas::cas<T, 7, 9>(data);
    avx2_cas::cas<T, 7, 8>(data);
    avx2_cas::cas<T, 5, 8>(data);
    avx2_cas::cas2<T>::template cas<5, 7, 6, 9>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas2<T>::template cas<6, 7, 0, 5>(data);
    avx2_cas::cas<T, 1, 6>(data);
    avx2_cas::cas4<T>::template cas<1, 5, 2, 7, 3, 8, 4, 9>(data);
    avx2_cas::cas2<T>::template cas<4, 8, 3, 7>(data);
    avx2_cas::cas3<T>::template cas<4, 7, 2, 5, 3, 6>(data);
    avx2_cas::cas2<T>::template cas<4, 6, 3, 5>(data);
    avx2_cas::cas3<T>::template cas<4, 5, 10, 11, 13, 14>(data);
    avx2_cas::cas<T, 12, 14>(data);
    avx2_cas::cas<T, 12, 13>(data);
    avx2_cas::cas<T, 10, 13>(data);
    avx2_cas::cas2<T>::template cas<10, 12, 11, 14>(data);
    avx2_cas::cas<T, 11, 13>(data);
    avx2_cas::cas2<T>::template cas<11, 12, 16, 17>(data);
    avx2_cas::cas<T, 15, 17>(data);
    avx2_cas::cas2<T>::template cas<15, 16, 19, 20>(data);
    avx2_cas::cas<T, 18, 20>(data);
    avx2_cas::cas<T, 18, 19>(data);
    avx2_cas::cas3<T>::template cas<15, 18, 16, 19, 17, 20>(data);
    avx2_cas::cas2<T>::template cas<17, 19, 16, 18>(data);
    avx2_cas::cas2<T>::template cas<17, 18, 10, 16>(data);
    avx2_cas::cas2<T>::template cas<10, 15, 11, 17>(data);
    avx2_cas::cas<T, 11, 16>(data);
    avx2_cas::cas4<T>::template cas<11, 15, 12, 18, 13, 19, 14, 20>(data);
    avx2_cas::cas2<T>::template cas<14, 19, 13, 18>(data);
    avx2_cas::cas3<T>::template cas<14, 18, 12, 15, 13, 16>(data);
    avx2_cas::cas<T, 14, 17>(data);
    avx2_cas::cas2<T>::template cas<14, 16, 13, 15>(data);
    avx2_cas::cas2<T>::template cas<14, 15, 0, 11>(data);
    avx2_cas::cas2<T>::template cas<0, 10, 1, 12>(data);
    avx2_cas::cas<T, 1, 11>(data);
    avx2_cas::cas4<T>::template cas<1, 10, 2, 13, 3, 14, 4, 15>(data);
    avx2_cas::cas2<T>::template cas<4, 14, 3, 13>(data);
    avx2_cas::cas3<T>::template cas<4, 13, 2, 10, 3, 11>(data);
    avx2_cas::cas<T, 4, 12>(data);
    avx2_cas::cas2<T>::template cas<4, 11, 3, 10>(data);
    avx2_cas::cas3<T>::template cas<4, 10, 5, 16, 6, 17>(data);
    avx2_cas::cas4<T>::template cas<6, 16, 7, 18, 8, 19, 9, 20>(data);
    avx2_cas::cas2<T>::template cas<9, 19, 8, 18>(data);
    avx2_cas::cas3<T>::template cas<9, 18, 7, 16, 8, 17>(data);
    avx2_cas::cas2<T>::template cas<9, 17, 8, 16>(data);
    avx2_cas::cas2<T>::template cas<9, 16, 5, 11>(data);
    avx2_cas::cas2<T>::template cas<5, 10, 6, 12>(data);
    avx2_cas::cas<T, 6, 11>(data);
    avx2_cas::cas4<T>::template cas<6, 10, 7, 13, 8, 14, 9, 15>(data);
    avx2_cas::cas2<T>::template cas<9, 14, 8, 13>(data);
    avx2_cas::cas3<T>::template cas<9, 13, 7, 10, 8, 11>(data);
    avx2_cas::cas<T, 9, 12>(data);
    avx2_cas::cas2<T>::template cas<9, 11, 8, 10>(data);
    avx2_cas::cas<T, 9, 10>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<22>::sort(T * const data)
{
    avx2_cas::cas2<T>::template cas<0, 1, 3, 4>(data);
    avx2_cas::cas<T, 2, 4>(data);
    avx2_cas::cas<T, 2, 3>(data);
    avx2_cas::cas<T, 0, 3>(data);
    avx2_cas::cas2<T>::template cas<0, 2, 1, 4>(data);
    avx2_cas::cas<T, 1, 3>(data);
    avx2_cas::cas2<T>::template cas<1, 2, 6, 7>(data);
    avx2_cas::cas<T, 5, 7>(data);
    avx2_cas::cas2<T>::template cas<5, 6, 9, 10>(data);
    avx2_cas::cas<T, 8, 10>(data);
    avx2_cas::cas<T, 8, 9>(data);
    avx2_cas::cas3<T>::template cas<5, 8, 6, 9, 7, 10>(data);
    avx2_cas::cas2<T>::template cas<7, 9, 6, 8>(data);
    avx2_cas::cas2<T>::template cas<7, 8, 0, 6>(data);
    avx2_cas::cas2<T>::template cas<0, 5, 1, 7>(data);
    avx2_cas::cas<T, 1, 6>(data);
    avx2_cas::cas4<T>::template cas<1, 5, 2, 8, 3, 9, 4, 10>(data);
    avx2_cas::cas2<T>::template cas<4, 9, 3, 8>(data);
    avx2_cas::cas3<T>::template cas<4, 8, 2, 5, 3, 6>(data);
    avx2_cas::cas<T, 4, 7>(data);
    avx2_cas::cas2<T>::template cas<4, 6, 3, 5>(data);
    avx2_cas::cas3<T>::template cas<4, 5, 11, 12, 14, 15>(data);
    avx2_cas::cas<T, 13, 15>(data);
    avx2_cas::cas<T, 13, 14>(data);
    avx2_cas::cas<T, 11, 14>(data);
    avx2_cas::cas2<T>::template cas<11, 13, 12, 15>(data);
    avx2_cas::cas<T, 12, 14>(data);
    avx2_cas::cas2<T>::template cas<12, 13, 17, 18>(data);
    avx2_cas::cas<T, 16, 18>(data);
    avx2_cas::cas2<T>::template cas<16, 17, 20, 21>(data);
    avx2_cas::cas<T, 19, 21>(data);
    avx2_cas::cas<T, 19, 20>(data);
    avx2_cas::cas3<T>::template cas<16, 19, 17, 20, 18, 21>(data);
    avx2_cas::cas2<T>::template cas<18, 20, 17, 19>(data);
    avx2_cas::cas2<T>::template cas<18, 19, 11, 17>(data);
    avx2_cas::cas2<T>::template cas<11, 16, 12, 18>(data);
    avx2_cas::cas<T, 12, 17>(data);
    avx2_cas::cas4<T>::template cas<12, 16, 13, 19, 14, 20, 15, 21>(data);
    avx2_cas::cas2<T>::template cas<15, 20, 14, 19>(data);
    avx2_cas::cas3<T>::template cas<15, 19, 13, 16, 14, 17>(data);
    avx2_cas::cas<T, 15, 18>(data);
    avx2_cas::cas2<T>::template cas<15, 17, 14, 16>(data);
    avx2_cas::cas3<T>::template cas<15, 16, 0, 11, 1, 12>(data);
    avx2_cas::cas4<T>::template cas<1, 11, 2, 13, 3, 14, 4, 15>(data);
    avx2_cas::cas2<T>::template cas<4, 14, 3, 13>(data);
    avx2_cas::cas3<T>::template cas<4, 13, 2, 11, 3, 12>(data);
    avx2_cas::cas2<T>::template cas<4, 12, 3, 11>(data);
    avx2_cas::cas4<T>::template cas<4, 11, 5, 16, 6, 17, 7, 18>(data);
    avx2_cas::cas2<T>::template cas<7, 17, 6, 16>(data);
    avx2_cas::cas4<T>::template cas<7, 16, 8, 19, 9, 20, 10, 21>(data);
    avx2_cas::cas2<T>::template cas<10, 20, 9, 19>(data);
    avx2_cas::cas3<T>::template cas<10, 19, 8, 16, 9, 17>(data);
    avx2_cas::cas<T, 10, 18>(data);
    avx2_cas::cas2<T>::template cas<10, 17, 9, 16>(data);
    avx2_cas::cas4<T>::template cas<10, 16, 5, 11, 6, 12, 7, 13>(data);
    avx2_cas::cas2<T>::template cas<7, 12, 6, 11>(data);
    avx2_cas::cas3<T>::template cas<7, 11, 8, 14, 9, 15>(data);
    avx2_cas::cas2<T>::template cas<10, 15, 9, 14>(data);
    avx2_cas::cas3<T>::template cas<10, 14, 8, 11, 9, 12>(data);
    avx2_cas::cas<T, 10, 13>(data);
    avx2_cas::cas2<T>::template cas<10, 12, 9, 11>(data);
    avx2_cas::cas<T, 10, 11>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<23>::sort(T * const data)
{
    avx2_cas::cas2<T>::template cas<0, 1, 3, 4>(data);
    avx2_cas::cas<T, 2, 4>(data);
    avx2_cas::cas<T, 2, 3>(data);
    avx2_cas::cas<T, 0, 3>(data);
    avx2_cas::cas2<T>::template cas<0, 2, 1, 4>(data);
    avx2_cas::cas<T, 1, 3>(data);
    avx2_cas::cas2<T>::template cas<1, 2, 6, 7>(data);
    avx2_cas::cas<T, 5, 7>(data);
    avx2_cas::cas2<T>::template cas<5, 6, 9, 10>(data);
    avx2_cas::cas<T, 8, 10>(data);
    avx2_cas::cas<T, 8, 9>(data);
    avx2_cas::cas3<T>::template cas<5, 8, 6, 9, 7, 10>(data);
    avx2_cas::cas2<T>::template cas<7, 9, 6, 8>(data);
    avx2_cas::cas2<T>::template cas<7, 8, 0, 6>(data);
    avx2_cas::cas2<T>::template cas<0, 5, 1, 7>(data);
    avx2_cas::cas<T, 1, 6>(data);
    avx2_cas::cas4<T>::template cas<1, 5, 2, 8, 3, 9, 4, 10>(data);
    avx2_cas::cas2<T>::template cas<4, 9, 3, 8>(data);
    avx2_cas::cas3<T>::template cas<4, 8, 2, 5, 3, 6>(data);
    avx2_cas::cas<T, 4, 7>(data);
    avx2_cas::cas2<T>::template cas<4, 6, 3, 5>(data);
    avx2_cas::cas2<T>::template cas<4, 5, 12, 13>(data);
    avx2_cas::cas<T, 11, 13>(data);
    avx2_cas::cas2<T>::template cas<11, 12, 15, 16>(data);
    avx2_cas::cas<T, 14, 16>(data);
    avx2_cas::cas<T, 14, 15>(data);
    avx2_cas::cas3<T>::template cas<11, 14, 12, 15, 13, 16>(data);
    avx2_cas::cas2<T>::template cas<13, 15, 12, 14>(data);
    avx2_cas::cas2<T>::template cas<13, 14, 18, 19>(data);
    avx2_cas::cas<T, 17, 19>(data);
    avx2_cas::cas2<T>::template cas<17, 18, 21, 22>(data);
    avx2_cas::cas<T, 20, 22>(data);
    avx2_cas::cas<T, 20, 21>(data);
    avx2_cas::cas3<T>::template cas<17, 20, 18, 21, 19, 22>(data);
    avx2_cas::cas2<T>::template cas<19, 21, 18, 20>(data);
    avx2_cas::cas3<T>::template cas<19, 20, 11, 17, 12, 18>(data);
    avx2_cas::cas<T, 13, 19>(data);
    avx2_cas::cas2<T>::template cas<13, 18, 12, 17>(data);
    avx2_cas::cas4<T>::template cas<13, 17, 14, 20, 15, 21, 16, 22>(data);
    avx2_cas::cas2<T>::template cas<16, 21, 15, 20>(data);
    avx2_cas::cas3<T>::template cas<16, 20, 14, 17, 15, 18>(data);
    avx2_cas::cas<T, 16, 19>(data);
    avx2_cas::cas2<T>::template cas<16, 18, 15, 17>(data);
    avx2_cas::cas2<T>::template cas<16, 17, 0, 12>(data);
    avx2_cas::cas2<T>::template cas<0, 11, 1, 13>(data);
    avx2_cas::cas<T, 1, 12>(data);
    avx2_cas::cas4<T>::template cas<1, 11, 2, 14, 3, 15, 4, 16>(data);
    avx2_cas::cas2<T>::template cas<4, 15, 3, 14>(data);
    avx2_cas::cas3<T>::template cas<4, 14, 2, 11, 3, 12>(data);
    avx2_cas::cas<T, 4, 13>(data);
    avx2_cas::cas2<T>::template cas<4, 12, 3, 11>(data);
    avx2_cas::cas4<T>::template cas<4, 11, 5, 17, 6, 18, 7, 19>(data);
    avx2_cas::cas2<T>::template cas<7, 18, 6, 17>(data);
    avx2_cas::cas4<T>::template cas<7, 17, 8, 20, 9, 21, 10, 22>(data);
    avx2_cas::cas2<T>::template cas<10, 21, 9, 20>(data);
    avx2_cas::cas3<T>::template cas<10, 20, 8, 17, 9, 18>(data);
    avx2_cas::cas<T, 10, 19>(data);
    avx2_cas::cas2<T>::template cas<10, 18, 9, 17>(data);
    avx2_cas::cas4<T>::template cas<10, 17, 5, 11, 6, 12, 7, 13>(data);
    avx2_cas::cas2<T>::template cas<7, 12, 6, 11>(data);
    avx2_cas::cas4<T>::template cas<7, 11, 8, 14, 9, 15, 10, 16>(data);
    avx2_cas::cas2<T>::template cas<10, 15, 9, 14>(data);
    avx2_cas::cas3<T>::template cas<10, 14, 8, 11, 9, 12>(data);
    avx2_cas::cas<T, 10, 13>(data);
    avx2_cas::cas2<T>::template cas<10, 12, 9, 11>(data);
    avx2_cas::cas<T, 10, 11>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<24>::sort(T * const data)
{
    avx2_cas::cas<T, 1, 2>(data);
    avx2_cas::cas<T, 0, 2>(data);
    avx2_cas::cas2<T>::template cas<0, 1, 4, 5>(data);
    avx2_cas::cas<T, 3, 5>(data);
    avx2_cas::cas<T, 3, 4>(data);
    avx2_cas::cas3<T>::template cas<0, 3, 1, 4, 2, 5>(data);
    avx2_cas::cas2<T>::template cas<2, 4, 1, 3>(data);
    avx2_cas::cas2<T>::template cas<2, 3, 7, 8>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas2<T>::template cas<6, 7, 10, 11>(data);
    avx2_cas::cas<T, 9, 11>(data);
    avx2_cas::cas<T, 9, 10>(data);
    avx2_cas::cas3<T>::template cas<6, 9, 7, 10, 8, 11>(data);
    avx2_cas::cas2<T>::template cas<8, 10, 7, 9>(data);
    avx2_cas::cas3<T>::template cas<8, 9, 0, 6, 1, 7>(data);
    avx2_cas::cas<T, 2, 8>(data);
    avx2_cas::cas2<T>::template cas<2, 7, 1, 6>(data);
    avx2_cas::cas4<T>::template cas<2, 6, 3, 9, 4, 10, 5, 11>(data);
    avx2_cas::cas2<T>::template cas<5, 10, 4, 9>(data);
    avx2_cas::cas3<T>::template cas<5, 9, 3, 6, 4, 7>(data);
    avx2_cas::cas<T, 5, 8>(data);
    avx2_cas::cas2<T>::template cas<5, 7, 4, 6>(data);
    avx2_cas::cas2<T>::template cas<5, 6, 13, 14>(data);
    avx2_cas::cas<T, 12, 14>(data);
    avx2_cas::cas2<T>::template cas<12, 13, 16, 17>(data);
    avx2_cas::cas<T, 15, 17>(data);
    avx2_cas::cas<T, 15, 16>(data);
    avx2_cas::cas3<T>::template cas<12, 15, 13, 16, 14, 17>(data);
    avx2_cas::cas2<T>::template cas<14, 16, 13, 15>(data);
    avx2_cas::cas2<T>::template cas<14, 15, 19, 20>(data);
    avx2_cas::cas<T, 18, 20>(data);
    avx2_cas::cas2<T>::template cas<18, 19, 22, 23>(data);
    avx2_cas::cas<T, 21, 23>(data);
    avx2_cas::cas<T, 21, 22>(data);
    avx2_cas::cas3<T>::template cas<18, 21, 19, 22, 20, 23>(data);
    avx2_cas::cas2<T>::template cas<20, 22, 19, 21>(data);
    avx2_cas::cas3<T>::template cas<20, 21, 12, 18, 13, 19>(data);
    avx2_cas::cas<T, 14, 20>(data);
    avx2_cas::cas2<T>::template cas<14, 19, 13, 18>(data);
    avx2_cas::cas4<T>::template cas<14, 18, 15, 21, 16, 22, 17, 23>(data);
    avx2_cas::cas2<T>::template cas<17, 22, 16, 21>(data);
    avx2_cas::cas3<T>::template cas<17, 21, 15, 18, 16, 19>(data);
    avx2_cas::cas<T, 17, 20>(data);
    avx2_cas::cas2<T>::template cas<17, 19, 16, 18>(data);
    avx2_cas::cas4<T>::template cas<17, 18, 0, 12, 1, 13, 2, 14>(data);
    avx2_cas::cas2<T>::template cas<2, 13, 1, 12>(data);
    avx2_cas::cas4<T>::template cas<2, 12, 3, 15, 4, 16, 5, 17>(data);
    avx2_cas::cas2<T>::template cas<5, 16, 4, 15>(data);
    avx2_cas::cas3<T>::template cas<5, 15, 3, 12, 4, 13>(data);
    avx2_cas::cas<T, 5, 14>(data);
    avx2_cas::cas2<T>::template cas<5, 13, 4, 12>(data);
    avx2_cas::cas4<T>::template cas<5, 12, 6, 18, 7, 19, 8, 20>(data);
    avx2_cas::cas2<T>::template cas<8, 19, 7, 18>(data);
    avx2_cas::cas4<T>::template cas<8, 18, 9, 21, 10, 22, 11, 23>(data);
    avx2_cas::cas2<T>::template cas<11, 22, 10, 21>(data);
    avx2_cas::cas3<T>::template cas<11, 21, 9, 18, 10, 19>(data);
    avx2_cas::cas<T, 11, 20>(data);
    avx2_cas::cas2<T>::template cas<11, 19, 10, 18>(data);
    avx2_cas::cas4<T>::template cas<11, 18, 6, 12, 7, 13, 8, 14>(data);
    avx2_cas::cas2<T>::template cas<8, 13, 7, 12>(data);
    avx2_cas::cas4<T>::template cas<8, 12, 9, 15, 10, 16, 11, 17>(data);
    avx2_cas::cas2<T>::template cas<11, 16, 10, 15>(data);
    avx2_cas::cas3<T>::template cas<11, 15, 9, 12, 10, 13>(data);
    avx2_cas::cas<T, 11, 14>(data);
    avx2_cas::cas2<T>::template cas<11, 13, 10, 12>(data);
    avx2_cas::cas<T, 11, 12>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<25>::sort(T * const data)
{
    avx2_cas::cas<T, 1, 2>(data);
    avx2_cas::cas<T, 0, 2>(data);
    avx2_cas::cas2<T>::template cas<0, 1, 4, 5>(data);
    avx2_cas::cas<T, 3, 5>(data);
    avx2_cas::cas<T, 3, 4>(data);
    avx2_cas::cas3<T>::template cas<0, 3, 1, 4, 2, 5>(data);
    avx2_cas::cas2<T>::template cas<2, 4, 1, 3>(data);
    avx2_cas::cas2<T>::template cas<2, 3, 7, 8>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas2<T>::template cas<6, 7, 10, 11>(data);
    avx2_cas::cas<T, 9, 11>(data);
    avx2_cas::cas<T, 9, 10>(data);
    avx2_cas::cas3<T>::template cas<6, 9, 7, 10, 8, 11>(data);
    avx2_cas::cas2<T>::template cas<8, 10, 7, 9>(data);
    avx2_cas::cas3<T>::template cas<8, 9, 0, 6, 1, 7>(data);
    avx2_cas::cas<T, 2, 8>(data);
    avx2_cas::cas2<T>::template cas<2, 7, 1, 6>(data);
    avx2_cas::cas4<T>::template cas<2, 6, 3, 9, 4, 10, 5, 11>(data);
    avx2_cas::cas2<T>::template cas<5, 10, 4, 9>(data);
    avx2_cas::cas3<T>::template cas<5, 9, 3, 6, 4, 7>(data);
    avx2_cas::cas<T, 5, 8>(data);
    avx2_cas::cas2<T>::template cas<5, 7, 4, 6>(data);
    avx2_cas::cas2<T>::template cas<5, 6, 13, 14>(data);
    avx2_cas::cas<T, 12, 14>(data);
    avx2_cas::cas2<T>::template cas<12, 13, 16, 17>(data);
    avx2_cas::cas<T, 15, 17>(data);
    avx2_cas::cas<T, 15, 16>(data);
    avx2_cas::cas3<T>::template cas<12, 15, 13, 16, 14, 17>(data);
    avx2_cas::cas2<T>::template cas<14, 16, 13, 15>(data);
    avx2_cas::cas2<T>::template cas<14, 15, 19, 20>(data);
    avx2_cas::cas<T, 18, 20>(data);
    avx2_cas::cas3<T>::template cas<18, 19, 21, 22, 23, 24>(data);
    avx2_cas::cas2<T>::template cas<21, 23, 22, 24>(data);
    avx2_cas::cas<T, 22, 23>(data);
    avx2_cas::cas<T, 18, 22>(data);
    avx2_cas::cas3<T>::template cas<18, 21, 19, 23, 20, 24>(data);
    avx2_cas::cas2<T>::template cas<20, 23, 19, 21>(data);
    avx2_cas::cas<T, 20, 22>(data);
    avx2_cas::cas2<T>::template cas<20, 21, 12, 19>(data);
    avx2_cas::cas3<T>::template cas<12, 18, 13, 20, 14, 21>(data);
    avx2_cas::cas2<T>::template cas<14, 20, 13, 18>(data);
    avx2_cas::cas<T, 14, 19>(data);
    avx2_cas::cas4<T>::template cas<14, 18, 15, 22, 16, 23, 17, 24>(data);
    avx2_cas::cas2<T>::template cas<17, 23, 16, 22>(data);
    avx2_cas::cas2<T>::template cas<17, 22, 15, 19>(data);
    avx2_cas::cas3<T>::template cas<15, 18, 16, 20, 17, 21>(data);
    avx2_cas::cas2<T>::template cas<17, 20, 16, 18>(data);
    avx2_cas::cas<T, 17, 19>(data);
    avx2_cas::cas2<T>::template cas<17, 18, 0, 13>(data);
    avx2_cas::cas3<T>::template cas<0, 12, 1, 14, 2, 15>(data);
    avx2_cas::cas2<T>::template cas<2, 14, 1, 12>(data);
    avx2_cas::cas<T, 2, 13>(data);
    avx2_cas::cas4<T>::template cas<2, 12, 3, 16, 4, 17, 5, 18>(data);
    avx2_cas::cas2<T>::template cas<5, 17, 4, 16>(data);
    avx2_cas::cas2<T>::template cas<5, 16, 3, 13>(data);
    avx2_cas::cas3<T>::template cas<3, 12, 4, 14, 5, 15>(data);
    avx2_cas::cas2<T>::template cas<5, 14, 4, 12>(data);
    avx2_cas::cas<T, 5, 13>(data);
    avx2_cas::cas4<T>::template cas<5, 12, 6, 19, 7, 20, 8, 21>(data);
    avx2_cas::cas2<T>::template cas<8, 20, 7, 19>(data);
    avx2_cas::cas4<T>::template cas<8, 19, 9, 22, 10, 23, 11, 24>(data);
    avx2_cas::cas2<T>::template cas<11, 23, 10, 22>(data);
    avx2_cas::cas3<T>::template cas<11, 22, 9, 19, 10, 20>(data);
    avx2_cas::cas<T, 11, 21>(data);
    avx2_cas::cas2<T>::template cas<11, 20, 10, 19>(data);
    avx2_cas::cas2<T>::template cas<11, 19, 6, 13>(data);
    avx2_cas::cas3<T>::template cas<6, 12, 7, 14, 8, 15>(data);
    avx2_cas::cas2<T>::template cas<8, 14, 7, 12>(data);
    avx2_cas::cas<T, 8, 13>(data);
    avx2_cas::cas4<T>::template cas<8, 12, 9, 16, 10, 17, 11, 18>(data);
    avx2_cas::cas2<T>::template cas<11, 17, 10, 16>(data);
    avx2_cas::cas2<T>::template cas<11, 16, 9, 13>(data);
    avx2_cas::cas3<T>::template cas<9, 12, 10, 14, 11, 15>(data);
    avx2_cas::cas2<T>::template cas<11, 14, 10, 12>(data);
    avx2_cas::cas<T, 11, 13>(data);
    avx2_cas::cas<T, 11, 12>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<26>::sort(T * const data)
{
    avx2_cas::cas<T, 1, 2>(data);
    avx2_cas::cas<T, 0, 2>(data);
    avx2_cas::cas2<T>::template cas<0, 1, 4, 5>(data);
    avx2_cas::cas<T, 3, 5>(data);
    avx2_cas::cas<T, 3, 4>(data);
    avx2_cas::cas3<T>::template cas<0, 3, 1, 4, 2, 5>(data);
    avx2_cas::cas2<T>::template cas<2, 4, 1, 3>(data);
    avx2_cas::cas2<T>::template cas<2, 3, 7, 8>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas3<T>::template cas<6, 7, 9, 10, 11, 12>(data);
    avx2_cas::cas2<T>::template cas<9, 11, 10, 12>(data);
    avx2_cas::cas<T, 10, 11>(data);
    avx2_cas::cas<T, 6, 10>(data);
    avx2_cas::cas3<T>::template cas<6, 9, 7, 11, 8, 12>(data);
    avx2_cas::cas2<T>::template cas<8, 11, 7, 9>(data);
    avx2_cas::cas<T, 8, 10>(data);
    avx2_cas::cas2<T>::template cas<8, 9, 0, 7>(data);
    avx2_cas::cas3<T>::template cas<0, 6, 1, 8, 2, 9>(data);
    avx2_cas::cas2<T>::template cas<2, 8, 1, 6>(data);
    avx2_cas::cas<T, 2, 7>(data);
    avx2_cas::cas4<T>::template cas<2, 6, 3, 10, 4, 11, 5, 12>(data);
    avx2_cas::cas2<T>::template cas<5, 11, 4, 10>(data);
    avx2_cas::cas2<T>::template cas<5, 10, 3, 7>(data);
    avx2_cas::cas3<T>::template cas<3, 6, 4, 8, 5, 9>(data);
    avx2_cas::cas2<T>::template cas<5, 8, 4, 6>(data);
    avx2_cas::cas<T, 5, 7>(data);
    avx2_cas::cas2<T>::template cas<5, 6, 14, 15>(data);
    avx2_cas::cas<T, 13, 15>(data);
    avx2_cas::cas2<T>::template cas<13, 14, 17, 18>(data);
    avx2_cas::cas<T, 16, 18>(data);
    avx2_cas::cas<T, 16, 17>(data);
    avx2_cas::cas3<T>::template cas<13, 16, 14, 17, 15, 18>(data);
    avx2_cas::cas2<T>::template cas<15, 17, 14, 16>(data);
    avx2_cas::cas2<T>::template cas<15, 16, 20, 21>(data);
    avx2_cas::cas<T, 19, 21>(data);
    avx2_cas::cas3<T>::template cas<19, 20, 22, 23, 24, 25>(data);
    avx2_cas::cas2<T>::template cas<22, 24, 23, 25>(data);
    avx2_cas::cas<T, 23, 24>(data);
    avx2_cas::cas<T, 19, 23>(data);
    avx2_cas::cas3<T>::template cas<19, 22, 20, 24, 21, 25>(data);
    avx2_cas::cas2<T>::template cas<21, 24, 20, 22>(data);
    avx2_cas::cas<T, 21, 23>(data);
    avx2_cas::cas2<T>::template cas<21, 22, 13, 20>(data);
    avx2_cas::cas3<T>::template cas<13, 19, 14, 21, 15, 22>(data);
    avx2_cas::cas2<T>::template cas<15, 21, 14, 19>(data);
    avx2_cas::cas<T, 15, 20>(data);
    avx2_cas::cas4<T>::template cas<15, 19, 16, 23, 17, 24, 18, 25>(data);
    avx2_cas::cas2<T>::template cas<18, 24, 17, 23>(data);
    avx2_cas::cas2<T>::template cas<18, 23, 16, 20>(data);
    avx2_cas::cas3<T>::template cas<16, 19, 17, 21, 18, 22>(data);
    avx2_cas::cas2<T>::template cas<18, 21, 17, 19>(data);
    avx2_cas::cas<T, 18, 20>(data);
    avx2_cas::cas4<T>::template cas<18, 19, 0, 13, 1, 14, 2, 15>(data);
    avx2_cas::cas2<T>::template cas<2, 14, 1, 13>(data);
    avx2_cas::cas4<T>::template cas<2, 13, 3, 16, 4, 17, 5, 18>(data);
    avx2_cas::cas2<T>::template cas<5, 17, 4, 16>(data);
    avx2_cas::cas3<T>::template cas<5, 16, 3, 13, 4, 14>(data);
    avx2_cas::cas<T, 5, 15>(data);
    avx2_cas::cas2<T>::template cas<5, 14, 4, 13>(data);
    avx2_cas::cas4<T>::template cas<5, 13, 6, 19, 7, 20, 8, 21>(data);
    avx2_cas::cas2<T>::template cas<8, 20, 7, 19>(data);
    avx2_cas::cas3<T>::template cas<8, 19, 9, 22, 10, 23>(data);
    avx2_cas::cas3<T>::template cas<10, 22, 11, 24, 12, 25>(data);
    avx2_cas::cas2<T>::template cas<12, 24, 11, 22>(data);
    avx2_cas::cas<T, 12, 23>(data);
    avx2_cas::cas3<T>::template cas<12, 22, 9, 19, 10, 20>(data);
    avx2_cas::cas2<T>::template cas<10, 19, 11, 21>(data);
    avx2_cas::cas2<T>::template cas<12, 21, 11, 19>(data);
    avx2_cas::cas<T, 12, 20>(data);
    avx2_cas::cas4<T>::template cas<12, 19, 6, 13, 7, 14, 8, 15>(data);
    avx2_cas::cas2<T>::template cas<8, 14, 7, 13>(data);
    avx2_cas::cas3<T>::template cas<8, 13, 9, 16, 10, 17>(data);
    avx2_cas::cas2<T>::template cas<10, 16, 11, 18>(data);
    avx2_cas::cas2<T>::template cas<12, 18, 11, 16>(data);
    avx2_cas::cas<T, 12, 17>(data);
    avx2_cas::cas3<T>::template cas<12, 16, 9, 13, 10, 14>(data);
    avx2_cas::cas2<T>::template cas<10, 13, 11, 15>(data);
    avx2_cas::cas2<T>::template cas<12, 15, 11, 13>(data);
    avx2_cas::cas<T, 12, 14>(data);
    avx2_cas::cas<T, 12, 13>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<27>::sort(T * const data)
{
    avx2_cas::cas<T, 1, 2>(data);
    avx2_cas::cas<T, 0, 2>(data);
    avx2_cas::cas2<T>::template cas<0, 1, 4, 5>(data);
    avx2_cas::cas<T, 3, 5>(data);
    avx2_cas::cas<T, 3, 4>(data);
    avx2_cas::cas3<T>::template cas<0, 3, 1, 4, 2, 5>(data);
    avx2_cas::cas2<T>::template cas<2, 4, 1, 3>(data);
    avx2_cas::cas2<T>::template cas<2, 3, 7, 8>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas3<T>::template cas<6, 7, 9, 10, 11, 12>(data);
    avx2_cas::cas2<T>::template cas<9, 11, 10, 12>(data);
    avx2_cas::cas<T, 10, 11>(data);
    avx2_cas::cas<T, 6, 10>(data);
    avx2_cas::cas3<T>::template cas<6, 9, 7, 11, 8, 12>(data);
    avx2_cas::cas2<T>::template cas<8, 11, 7, 9>(data);
    avx2_cas::cas<T, 8, 10>(data);
    avx2_cas::cas2<T>::template cas<8, 9, 0, 7>(data);
    avx2_cas::cas3<T>::template cas<0, 6, 1, 8, 2, 9>(data);
    avx2_cas::cas2<T>::template cas<2, 8, 1, 6>(data);
    avx2_cas::cas<T, 2, 7>(data);
    avx2_cas::cas4<T>::template cas<2, 6, 3, 10, 4, 11, 5, 12>(data);
    avx2_cas::cas2<T>::template cas<5, 11, 4, 10>(data);
    avx2_cas::cas2<T>::template cas<5, 10, 3, 7>(data);
    avx2_cas::cas3<T>::template cas<3, 6, 4, 8, 5, 9>(data);
    avx2_cas::cas2<T>::template cas<5, 8, 4, 6>(data);
    avx2_cas::cas<T, 5, 7>(data);
    avx2_cas::cas2<T>::template cas<5, 6, 14, 15>(data);
    avx2_cas::cas<T, 13, 15>(data);
    avx2_cas::cas3<T>::template cas<13, 14, 16, 17, 18, 19>(data);
    avx2_cas::cas2<T>::template cas<16, 18, 17, 19>(data);
    avx2_cas::cas<T, 17, 18>(data);
    avx2_cas::cas<T, 13, 17>(data);
    avx2_cas::cas3<T>::template cas<13, 16, 14, 18, 15, 19>(data);
    avx2_cas::cas2<T>::template cas<15, 18, 14, 16>(data);
    avx2_cas::cas<T, 15, 17>(data);
    avx2_cas::cas2<T>::template cas<15, 16, 21, 22>(data);
    avx2_cas::cas<T, 20, 22>(data);
    avx2_cas::cas3<T>::template cas<20, 21, 23, 24, 25, 26>(data);
    avx2_cas::cas2<T>::template cas<23, 25, 24, 26>(data);
    avx2_cas::cas<T, 24, 25>(data);
    avx2_cas::cas<T, 20, 24>(data);
    avx2_cas::cas3<T>::template cas<20, 23, 21, 25, 22, 26>(data);
    avx2_cas::cas2<T>::template cas<22, 25, 21, 23>(data);
    avx2_cas::cas<T, 22, 24>(data);
    avx2_cas::cas3<T>::template cas<22, 23, 13, 20, 14, 21>(data);
    avx2_cas::cas<T, 15, 22>(data);
    avx2_cas::cas2<T>::template cas<15, 21, 14, 20>(data);
    avx2_cas::cas3<T>::template cas<15, 20, 16, 23, 17, 24>(data);
    avx2_cas::cas3<T>::template cas<17, 23, 18, 25, 19, 26>(data);
    avx2_cas::cas2<T>::template cas<19, 25, 18, 23>(data);
    avx2_cas::cas<T, 19, 24>(data);
    avx2_cas::cas3<T>::template cas<19, 23, 16, 20, 17, 21>(data);
    avx2_cas::cas2<T>::template cas<17, 20, 18, 22>(data);
    avx2_cas::cas2<T>::template cas<19, 22, 18, 20>(data);
    avx2_cas::cas<T, 19, 21>(data);
    avx2_cas::cas2<T>::template cas<19, 20, 0, 14>(data);
    avx2_cas::cas3<T>::template cas<0, 13, 1, 15, 2, 16>(data);
    avx2_cas::cas2<T>::template cas<2, 15, 1, 13>(data);
    avx2_cas::cas<T, 2, 14>(data);
    avx2_cas::cas4<T>::template cas<2, 13, 3, 17, 4, 18, 5, 19>(data);
    avx2_cas::cas2<T>::template cas<5, 18, 4, 17>(data);
    avx2_cas::cas2<T>::template cas<5, 17, 3, 14>(data);
    avx2_cas::cas3<T>::template cas<3, 13, 4, 15, 5, 16>(data);
    avx2_cas::cas2<T>::template cas<5, 15, 4, 13>(data);
    avx2_cas::cas<T, 5, 14>(data);
    avx2_cas::cas4<T>::template cas<5, 13, 6, 20, 7, 21, 8, 22>(data);
    avx2_cas::cas2<T>::template cas<8, 21, 7, 20>(data);
    avx2_cas::cas3<T>::template cas<8, 20, 9, 23, 10, 24>(data);
    avx2_cas::cas3<T>::template cas<10, 23, 11, 25, 12, 26>(data);
    avx2_cas::cas2<T>::template cas<12, 25, 11, 23>(data);
    avx2_cas::cas<T, 12, 24>(data);
    avx2_cas::cas3<T>::template cas<12, 23, 9, 20, 10, 21>(data);
    avx2_cas::cas2<T>::template cas<10, 20, 11, 22>(data);
    avx2_cas::cas2<T>::template cas<12, 22, 11, 20>(data);
    avx2_cas::cas<T, 12, 21>(data);
    avx2_cas::cas4<T>::template cas<12, 20, 6, 13, 7, 14, 8, 15>(data);
    avx2_cas::cas2<T>::template cas<8, 14, 7, 13>(data);
    avx2_cas::cas3<T>::template cas<8, 13, 9, 16, 10, 17>(data);
    avx2_cas::cas3<T>::template cas<10, 16, 11, 18, 12, 19>(data);
    avx2_cas::cas2<T>::template cas<12, 18, 11, 16>(data);
    avx2_cas::cas<T, 12, 17>(data);
    avx2_cas::cas3<T>::template cas<12, 16, 9, 13, 10, 14>(data);
    avx2_cas::cas2<T>::template cas<10, 13, 11, 15>(data);
    avx2_cas::cas2<T>::template cas<12, 15, 11, 13>(data);
    avx2_cas::cas<T, 12, 14>(data);
    avx2_cas::cas<T, 12, 13>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<28>::sort(T * const data)
{
    avx2_cas::cas<T, 1, 2>(data);
    avx2_cas::cas<T, 0, 2>(data);
    avx2_cas::cas3<T>::template cas<0, 1, 3, 4, 5, 6>(data);
    avx2_cas::cas2<T>::template cas<3, 5, 4, 6>(data);
    avx2_cas::cas<T, 4, 5>(data);
    avx2_cas::cas<T, 0, 4>(data);
    avx2_cas::cas3<T>::template cas<0, 3, 1, 5, 2, 6>(data);
    avx2_cas::cas2<T>::template cas<2, 5, 1, 3>(data);
    avx2_cas::cas<T, 2, 4>(data);
    avx2_cas::cas2<T>::template cas<2, 3, 8, 9>(data);
    avx2_cas::cas<T, 7, 9>(data);
    avx2_cas::cas3<T>::template cas<7, 8, 10, 11, 12, 13>(data);
    avx2_cas::cas2<T>::template cas<10, 12, 11, 13>(data);
    avx2_cas::cas<T, 11, 12>(data);
    avx2_cas::cas<T, 7, 11>(data);
    avx2_cas::cas3<T>::template cas<7, 10, 8, 12, 9, 13>(data);
    avx2_cas::cas2<T>::template cas<9, 12, 8, 10>(data);
    avx2_cas::cas<T, 9, 11>(data);
    avx2_cas::cas3<T>::template cas<9, 10, 0, 7, 1, 8>(data);
    avx2_cas::cas<T, 2, 9>(data);
    avx2_cas::cas2<T>::template cas<2, 8, 1, 7>(data);
    avx2_cas::cas3<T>::template cas<2, 7, 3, 10, 4, 11>(data);
    avx2_cas::cas3<T>::template cas<4, 10, 5, 12, 6, 13>(data);
    avx2_cas::cas2<T>::template cas<6, 12, 5, 10>(data);
    avx2_cas::cas<T, 6, 11>(data);
    avx2_cas::cas3<T>::template cas<6, 10, 3, 7, 4, 8>(data);
    avx2_cas::cas2<T>::template cas<4, 7, 5, 9>(data);
    avx2_cas::cas2<T>::template cas<6, 9, 5, 7>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas2<T>::template cas<6, 7, 15, 16>(data);
    avx2_cas::cas<T, 14, 16>(data);
    avx2_cas::cas3<T>::template cas<14, 15, 17, 18, 19, 20>(data);
    avx2_cas::cas2<T>::template cas<17, 19, 18, 20>(data);
    avx2_cas::cas<T, 18, 19>(data);
    avx2_cas::cas<T, 14, 18>(data);
    avx2_cas::cas3<T>::template cas<14, 17, 15, 19, 16, 20>(data);
    avx2_cas::cas2<T>::template cas<16, 19, 15, 17>(data);
    avx2_cas::cas<T, 16, 18>(data);
    avx2_cas::cas2<T>::template cas<16, 17, 22, 23>(data);
    avx2_cas::cas<T, 21, 23>(data);
    avx2_cas::cas3<T>::template cas<21, 22, 24, 25, 26, 27>(data);
    avx2_cas::cas2<T>::template cas<24, 26, 25, 27>(data);
    avx2_cas::cas<T, 25, 26>(data);
    avx2_cas::cas<T, 21, 25>(data);
    avx2_cas::cas3<T>::template cas<21, 24, 22, 26, 23, 27>(data);
    avx2_cas::cas2<T>::template cas<23, 26, 22, 24>(data);
    avx2_cas::cas<T, 23, 25>(data);
    avx2_cas::cas3<T>::template cas<23, 24, 14, 21, 15, 22>(data);
    avx2_cas::cas<T, 16, 23>(data);
    avx2_cas::cas2<T>::template cas<16, 22, 15, 21>(data);
    avx2_cas::cas3<T>::template cas<16, 21, 17, 24, 18, 25>(data);
    avx2_cas::cas3<T>::template cas<18, 24, 19, 26, 20, 27>(data);
    avx2_cas::cas2<T>::template cas<20, 26, 19, 24>(data);
    avx2_cas::cas<T, 20, 25>(data);
    avx2_cas::cas3<T>::template cas<20, 24, 17, 21, 18, 22>(data);
    avx2_cas::cas2<T>::template cas<18, 21, 19, 23>(data);
    avx2_cas::cas2<T>::template cas<20, 23, 19, 21>(data);
    avx2_cas::cas<T, 20, 22>(data);
    avx2_cas::cas4<T>::template cas<20, 21, 0, 14, 1, 15, 2, 16>(data);
    avx2_cas::cas2<T>::template cas<2, 15, 1, 14>(data);
    avx2_cas::cas3<T>::template cas<2, 14, 3, 17, 4, 18>(data);
    avx2_cas::cas3<T>::template cas<4, 17, 5, 19, 6, 20>(data);
    avx2_cas::cas2<T>::template cas<6, 19, 5, 17>(data);
    avx2_cas::cas<T, 6, 18>(data);
    avx2_cas::cas3<T>::template cas<6, 17, 3, 14, 4, 15>(data);
    avx2_cas::cas2<T>::template cas<4, 14, 5, 16>(data);
    avx2_cas::cas2<T>::template cas<6, 16, 5, 14>(data);
    avx2_cas::cas<T, 6, 15>(data);
    avx2_cas::cas4<T>::template cas<6, 14, 7, 21, 8, 22, 9, 23>(data);
    avx2_cas::cas2<T>::template cas<9, 22, 8, 21>(data);
    avx2_cas::cas3<T>::template cas<9, 21, 10, 24, 11, 25>(data);
    avx2_cas::cas3<T>::template cas<11, 24, 12, 26, 13, 27>(data);
    avx2_cas::cas2<T>::template cas<13, 26, 12, 24>(data);
    avx2_cas::cas<T, 13, 25>(data);
    avx2_cas::cas3<T>::template cas<13, 24, 10, 21, 11, 22>(data);
    avx2_cas::cas2<T>::template cas<11, 21, 12, 23>(data);
    avx2_cas::cas2<T>::template cas<13, 23, 12, 21>(data);
    avx2_cas::cas<T, 13, 22>(data);
    avx2_cas::cas4<T>::template cas<13, 21, 7, 14, 8, 15, 9, 16>(data);
    avx2_cas::cas2<T>::template cas<9, 15, 8, 14>(data);
    avx2_cas::cas3<T>::template cas<9, 14, 10, 17, 11, 18>(data);
    avx2_cas::cas3<T>::template cas<11, 17, 12, 19, 13, 20>(data);
    avx2_cas::cas2<T>::template cas<13, 19, 12, 17>(data);
    avx2_cas::cas<T, 13, 18>(data);
    avx2_cas::cas3<T>::template cas<13, 17, 10, 14, 11, 15>(data);
    avx2_cas::cas2<T>::template cas<11, 14, 12, 16>(data);
    avx2_cas::cas2<T>::template cas<13, 16, 12, 14>(data);
    avx2_cas::cas<T, 13, 15>(data);
    avx2_cas::cas<T, 13, 14>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<29>::sort(T * const data)
{
    avx2_cas::cas<T, 1, 2>(data);
    avx2_cas::cas<T, 0, 2>(data);
    avx2_cas::cas3<T>::template cas<0, 1, 3, 4, 5, 6>(data);
    avx2_cas::cas2<T>::template cas<3, 5, 4, 6>(data);
    avx2_cas::cas<T, 4, 5>(data);
    avx2_cas::cas<T, 0, 4>(data);
    avx2_cas::cas3<T>::template cas<0, 3, 1, 5, 2, 6>(data);
    avx2_cas::cas2<T>::template cas<2, 5, 1, 3>(data);
    avx2_cas::cas<T, 2, 4>(data);
    avx2_cas::cas2<T>::template cas<2, 3, 8, 9>(data);
    avx2_cas::cas<T, 7, 9>(data);
    avx2_cas::cas3<T>::template cas<7, 8, 10, 11, 12, 13>(data);
    avx2_cas::cas2<T>::template cas<10, 12, 11, 13>(data);
    avx2_cas::cas<T, 11, 12>(data);
    avx2_cas::cas<T, 7, 11>(data);
    avx2_cas::cas3<T>::template cas<7, 10, 8, 12, 9, 13>(data);
    avx2_cas::cas2<T>::template cas<9, 12, 8, 10>(data);
    avx2_cas::cas<T, 9, 11>(data);
    avx2_cas::cas3<T>::template cas<9, 10, 0, 7, 1, 8>(data);
    avx2_cas::cas<T, 2, 9>(data);
    avx2_cas::cas2<T>::template cas<2, 8, 1, 7>(data);
    avx2_cas::cas3<T>::template cas<2, 7, 3, 10, 4, 11>(data);
    avx2_cas::cas3<T>::template cas<4, 10, 5, 12, 6, 13>(data);
    avx2_cas::cas2<T>::template cas<6, 12, 5, 10>(data);
    avx2_cas::cas<T, 6, 11>(data);
    avx2_cas::cas3<T>::template cas<6, 10, 3, 7, 4, 8>(data);
    avx2_cas::cas2<T>::template cas<4, 7, 5, 9>(data);
    avx2_cas::cas2<T>::template cas<6, 9, 5, 7>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas2<T>::template cas<6, 7, 15, 16>(data);
    avx2_cas::cas<T, 14, 16>(data);
    avx2_cas::cas3<T>::template cas<14, 15, 17, 18, 19, 20>(data);
    avx2_cas::cas2<T>::template cas<17, 19, 18, 20>(data);
    avx2_cas::cas<T, 18, 19>(data);
    avx2_cas::cas<T, 14, 18>(data);
    avx2_cas::cas3<T>::template cas<14, 17, 15, 19, 16, 20>(data);
    avx2_cas::cas2<T>::template cas<16, 19, 15, 17>(data);
    avx2_cas::cas<T, 16, 18>(data);
    avx2_cas::cas3<T>::template cas<16, 17, 21, 22, 23, 24>(data);
    avx2_cas::cas2<T>::template cas<21, 23, 22, 24>(data);
    avx2_cas::cas3<T>::template cas<22, 23, 25, 26, 27, 28>(data);
    avx2_cas::cas2<T>::template cas<25, 27, 26, 28>(data);
    avx2_cas::cas2<T>::template cas<26, 27, 21, 25>(data);
    avx2_cas::cas<T, 22, 26>(data);
    avx2_cas::cas3<T>::template cas<22, 25, 23, 27, 24, 28>(data);
    avx2_cas::cas2<T>::template cas<24, 27, 23, 25>(data);
    avx2_cas::cas<T, 24, 26>(data);
    avx2_cas::cas2<T>::template cas<24, 25, 14, 22>(data);
    avx2_cas::cas3<T>::template cas<14, 21, 15, 23, 16, 24>(data);
    avx2_cas::cas2<T>::template cas<16, 23, 15, 21>(data);
    avx2_cas::cas<T, 16, 22>(data);
    avx2_cas::cas3<T>::template cas<16, 21, 17, 25, 18, 26>(data);
    avx2_cas::cas3<T>::template cas<18, 25, 19, 27, 20, 28>(data);
    avx2_cas::cas2<T>::template cas<20, 27, 19, 25>(data);
    avx2_cas::cas<T, 20, 26>(data);
    avx2_cas::cas3<T>::template cas<20, 25, 17, 21, 18, 22>(data);
    avx2_cas::cas3<T>::template cas<18, 21, 19, 23, 20, 24>(data);
    avx2_cas::cas2<T>::template cas<20, 23, 19, 21>(data);
    avx2_cas::cas<T, 20, 22>(data);
    avx2_cas::cas2<T>::template cas<20, 21, 0, 15>(data);
    avx2_cas::cas3<T>::template cas<0, 14, 1, 16, 2, 17>(data);
    avx2_cas::cas2<T>::template cas<2, 16, 1, 14>(data);
    avx2_cas::cas<T, 2, 15>(data);
    avx2_cas::cas3<T>::template cas<2, 14, 3, 18, 4, 19>(data);
    avx2_cas::cas3<T>::template cas<4, 18, 5, 20, 6, 21>(data);
    avx2_cas::cas2<T>::template cas<6, 20, 5, 18>(data);
    avx2_cas::cas<T, 6, 19>(data);
    avx2_cas::cas3<T>::template cas<6, 18, 3, 14, 4, 15>(data);
    avx2_cas::cas3<T>::template cas<4, 14, 5, 16, 6, 17>(data);
    avx2_cas::cas2<T>::template cas<6, 16, 5, 14>(data);
    avx2_cas::cas<T, 6, 15>(data);
    avx2_cas::cas4<T>::template cas<6, 14, 7, 22, 8, 23, 9, 24>(data);
    avx2_cas::cas2<T>::template cas<9, 23, 8, 22>(data);
    avx2_cas::cas3<T>::template cas<9, 22, 10, 25, 11, 26>(data);
    avx2_cas::cas3<T>::template cas<11, 25, 12, 27, 13, 28>(data);
    avx2_cas::cas2<T>::template cas<13, 27, 12, 25>(data);
    avx2_cas::cas<T, 13, 26>(data);
    avx2_cas::cas3<T>::template cas<13, 25, 10, 22, 11, 23>(data);
    avx2_cas::cas2<T>::template cas<11, 22, 12, 24>(data);
    avx2_cas::cas2<T>::template cas<13, 24, 12, 22>(data);
    avx2_cas::cas<T, 13, 23>(data);
    avx2_cas::cas2<T>::template cas<13, 22, 7, 15>(data);
    avx2_cas::cas3<T>::template cas<7, 14, 8, 16, 9, 17>(data);
    avx2_cas::cas2<T>::template cas<9, 16, 8, 14>(data);
    avx2_cas::cas<T, 9, 15>(data);
    avx2_cas::cas3<T>::template cas<9, 14, 10, 18, 11, 19>(data);
    avx2_cas::cas3<T>::template cas<11, 18, 12, 20, 13, 21>(data);
    avx2_cas::cas2<T>::template cas<13, 20, 12, 18>(data);
    avx2_cas::cas<T, 13, 19>(data);
    avx2_cas::cas3<T>::template cas<13, 18, 10, 14, 11, 15>(data);
    avx2_cas::cas3<T>::template cas<11, 14, 12, 16, 13, 17>(data);
    avx2_cas::cas2<T>::template cas<13, 16, 12, 14>(data);
    avx2_cas::cas<T, 13, 15>(data);
    avx2_cas::cas<T, 13, 14>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<30>::sort(T * const data)
{
    avx2_cas::cas<T, 1, 2>(data);
    avx2_cas::cas<T, 0, 2>(data);
    avx2_cas::cas3<T>::template cas<0, 1, 3, 4, 5, 6>(data);
    avx2_cas::cas2<T>::template cas<3, 5, 4, 6>(data);
    avx2_cas::cas<T, 4, 5>(data);
    avx2_cas::cas<T, 0, 4>(data);
    avx2_cas::cas3<T>::template cas<0, 3, 1, 5, 2, 6>(data);
    avx2_cas::cas2<T>::template cas<2, 5, 1, 3>(data);
    avx2_cas::cas<T, 2, 4>(data);
    avx2_cas::cas3<T>::template cas<2, 3, 7, 8, 9, 10>(data);
    avx2_cas::cas2<T>::template cas<7, 9, 8, 10>(data);
    avx2_cas::cas3<T>::template cas<8, 9, 11, 12, 13, 14>(data);
    avx2_cas::cas2<T>::template cas<11, 13, 12, 14>(data);
    avx2_cas::cas2<T>::template cas<12, 13, 7, 11>(data);
    avx2_cas::cas<T, 8, 12>(data);
    avx2_cas::cas3<T>::template cas<8, 11, 9, 13, 10, 14>(data);
    avx2_cas::cas2<T>::template cas<10, 13, 9, 11>(data);
    avx2_cas::cas<T, 10, 12>(data);
    avx2_cas::cas2<T>::template cas<10, 11, 0, 8>(data);
    avx2_cas::cas3<T>::template cas<0, 7, 1, 9, 2, 10>(data);
    avx2_cas::cas2<T>::template cas<2, 9, 1, 7>(data);
    avx2_cas::cas<T, 2, 8>(data);
    avx2_cas::cas3<T>::template cas<2, 7, 3, 11, 4, 12>(data);
    avx2_cas::cas3<T>::template cas<4, 11, 5, 13, 6, 14>(data);
    avx2_cas::cas2<T>::template cas<6, 13, 5, 11>(data);
    avx2_cas::cas<T, 6, 12>(data);
    avx2_cas::cas3<T>::template cas<6, 11, 3, 7, 4, 8>(data);
    avx2_cas::cas3<T>::template cas<4, 7, 5, 9, 6, 10>(data);
    avx2_cas::cas2<T>::template cas<6, 9, 5, 7>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas2<T>::template cas<6, 7, 16, 17>(data);
    avx2_cas::cas<T, 15, 17>(data);
    avx2_cas::cas3<T>::template cas<15, 16, 18, 19, 20, 21>(data);
    avx2_cas::cas2<T>::template cas<18, 20, 19, 21>(data);
    avx2_cas::cas<T, 19, 20>(data);
    avx2_cas::cas<T, 15, 19>(data);
    avx2_cas::cas3<T>::template cas<15, 18, 16, 20, 17, 21>(data);
    avx2_cas::cas2<T>::template cas<17, 20, 16, 18>(data);
    avx2_cas::cas<T, 17, 19>(data);
    avx2_cas::cas3<T>::template cas<17, 18, 22, 23, 24, 25>(data);
    avx2_cas::cas2<T>::template cas<22, 24, 23, 25>(data);
    avx2_cas::cas3<T>::template cas<23, 24, 26, 27, 28, 29>(data);
    avx2_cas::cas2<T>::template cas<26, 28, 27, 29>(data);
    avx2_cas::cas2<T>::template cas<27, 28, 22, 26>(data);
    avx2_cas::cas<T, 23, 27>(data);
    avx2_cas::cas3<T>::template cas<23, 26, 24, 28, 25, 29>(data);
    avx2_cas::cas2<T>::template cas<25, 28, 24, 26>(data);
    avx2_cas::cas<T, 25, 27>(data);
    avx2_cas::cas2<T>::template cas<25, 26, 15, 23>(data);
    avx2_cas::cas3<T>::template cas<15, 22, 16, 24, 17, 25>(data);
    avx2_cas::cas2<T>::template cas<17, 24, 16, 22>(data);
    avx2_cas::cas<T, 17, 23>(data);
    avx2_cas::cas3<T>::template cas<17, 22, 18, 26, 19, 27>(data);
    avx2_cas::cas3<T>::template cas<19, 26, 20, 28, 21, 29>(data);
    avx2_cas::cas2<T>::template cas<21, 28, 20, 26>(data);
    avx2_cas::cas<T, 21, 27>(data);
    avx2_cas::cas3<T>::template cas<21, 26, 18, 22, 19, 23>(data);
    avx2_cas::cas3<T>::template cas<19, 22, 20, 24, 21, 25>(data);
    avx2_cas::cas2<T>::template cas<21, 24, 20, 22>(data);
    avx2_cas::cas<T, 21, 23>(data);
    avx2_cas::cas4<T>::template cas<21, 22, 0, 15, 1, 16, 2, 17>(data);
    avx2_cas::cas2<T>::template cas<2, 16, 1, 15>(data);
    avx2_cas::cas3<T>::template cas<2, 15, 3, 18, 4, 19>(data);
    avx2_cas::cas3<T>::template cas<4, 18, 5, 20, 6, 21>(data);
    avx2_cas::cas2<T>::template cas<6, 20, 5, 18>(data);
    avx2_cas::cas<T, 6, 19>(data);
    avx2_cas::cas3<T>::template cas<6, 18, 3, 15, 4, 16>(data);
    avx2_cas::cas2<T>::template cas<4, 15, 5, 17>(data);
    avx2_cas::cas2<T>::template cas<6, 17, 5, 15>(data);
    avx2_cas::cas<T, 6, 16>(data);
    avx2_cas::cas3<T>::template cas<6, 15, 7, 22, 8, 23>(data);
    avx2_cas::cas3<T>::template cas<8, 22, 9, 24, 10, 25>(data);
    avx2_cas::cas2<T>::template cas<10, 24, 9, 22>(data);
    avx2_cas::cas<T, 10, 23>(data);
    avx2_cas::cas3<T>::template cas<10, 22, 11, 26, 12, 27>(data);
    avx2_cas::cas3<T>::template cas<12, 26, 13, 28, 14, 29>(data);
    avx2_cas::cas2<T>::template cas<14, 28, 13, 26>(data);
    avx2_cas::cas<T, 14, 27>(data);
    avx2_cas::cas3<T>::template cas<14, 26, 11, 22, 12, 23>(data);
    avx2_cas::cas3<T>::template cas<12, 22, 13, 24, 14, 25>(data);
    avx2_cas::cas2<T>::template cas<14, 24, 13, 22>(data);
    avx2_cas::cas<T, 14, 23>(data);
    avx2_cas::cas3<T>::template cas<14, 22, 7, 15, 8, 16>(data);
    avx2_cas::cas3<T>::template cas<8, 15, 9, 17, 10, 18>(data);
    avx2_cas::cas2<T>::template cas<10, 17, 9, 15>(data);
    avx2_cas::cas<T, 10, 16>(data);
    avx2_cas::cas3<T>::template cas<10, 15, 11, 19, 12, 20>(data);
    avx2_cas::cas2<T>::template cas<12, 19, 13, 21>(data);
    avx2_cas::cas2<T>::template cas<14, 21, 13, 19>(data);
    avx2_cas::cas<T, 14, 20>(data);
    avx2_cas::cas3<T>::template cas<14, 19, 11, 15, 12, 16>(data);
    avx2_cas::cas3<T>::template cas<12, 15, 13, 17, 14, 18>(data);
    avx2_cas::cas2<T>::template cas<14, 17, 13, 15>(data);
    avx2_cas::cas<T, 14, 16>(data);
    avx2_cas::cas<T, 14, 15>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<31>::sort(T * const data)
{
    avx2_cas::cas<T, 1, 2>(data);
    avx2_cas::cas<T, 0, 2>(data);
    avx2_cas::cas3<T>::template cas<0, 1, 3, 4, 5, 6>(data);
    avx2_cas::cas2<T>::template cas<3, 5, 4, 6>(data);
    avx2_cas::cas<T, 4, 5>(data);
    avx2_cas::cas<T, 0, 4>(data);
    avx2_cas::cas3<T>::template cas<0, 3, 1, 5, 2, 6>(data);
    avx2_cas::cas2<T>::template cas<2, 5, 1, 3>(data);
    avx2_cas::cas<T, 2, 4>(data);
    avx2_cas::cas3<T>::template cas<2, 3, 7, 8, 9, 10>(data);
    avx2_cas::cas2<T>::template cas<7, 9, 8, 10>(data);
    avx2_cas::cas3<T>::template cas<8, 9, 11, 12, 13, 14>(data);
    avx2_cas::cas2<T>::template cas<11, 13, 12, 14>(data);
    avx2_cas::cas2<T>::template cas<12, 13, 7, 11>(data);
    avx2_cas::cas<T, 8, 12>(data);
    avx2_cas::cas3<T>::template cas<8, 11, 9, 13, 10, 14>(data);
    avx2_cas::cas2<T>::template cas<10, 13, 9, 11>(data);
    avx2_cas::cas<T, 10, 12>(data);
    avx2_cas::cas2<T>::template cas<10, 11, 0, 8>(data);
    avx2_cas::cas3<T>::template cas<0, 7, 1, 9, 2, 10>(data);
    avx2_cas::cas2<T>::template cas<2, 9, 1, 7>(data);
    avx2_cas::cas<T, 2, 8>(data);
    avx2_cas::cas3<T>::template cas<2, 7, 3, 11, 4, 12>(data);
    avx2_cas::cas3<T>::template cas<4, 11, 5, 13, 6, 14>(data);
    avx2_cas::cas2<T>::template cas<6, 13, 5, 11>(data);
    avx2_cas::cas<T, 6, 12>(data);
    avx2_cas::cas3<T>::template cas<6, 11, 3, 7, 4, 8>(data);
    avx2_cas::cas3<T>::template cas<4, 7, 5, 9, 6, 10>(data);
    avx2_cas::cas2<T>::template cas<6, 9, 5, 7>(data);
    avx2_cas::cas<T, 6, 8>(data);
    avx2_cas::cas3<T>::template cas<6, 7, 15, 16, 17, 18>(data);
    avx2_cas::cas2<T>::template cas<15, 17, 16, 18>(data);
    avx2_cas::cas3<T>::template cas<16, 17, 19, 20, 21, 22>(data);
    avx2_cas::cas2<T>::template cas<19, 21, 20, 22>(data);
    avx2_cas::cas2<T>::template cas<20, 21, 15, 19>(data);
    avx2_cas::cas<T, 16, 20>(data);
    avx2_cas::cas3<T>::template cas<16, 19, 17, 21, 18, 22>(data);
    avx2_cas::cas2<T>::template cas<18, 21, 17, 19>(data);
    avx2_cas::cas<T, 18, 20>(data);
    avx2_cas::cas3<T>::template cas<18, 19, 23, 24, 25, 26>(data);
    avx2_cas::cas2<T>::template cas<23, 25, 24, 26>(data);
    avx2_cas::cas3<T>::template cas<24, 25, 27, 28, 29, 30>(data);
    avx2_cas::cas2<T>::template cas<27, 29, 28, 30>(data);
    avx2_cas::cas2<T>::template cas<28, 29, 23, 27>(data);
    avx2_cas::cas<T, 24, 28>(data);
    avx2_cas::cas3<T>::template cas<24, 27, 25, 29, 26, 30>(data);
    avx2_cas::cas2<T>::template cas<26, 29, 25, 27>(data);
    avx2_cas::cas<T, 26, 28>(data);
    avx2_cas::cas3<T>::template cas<26, 27, 15, 23, 16, 24>(data);
    avx2_cas::cas3<T>::template cas<16, 23, 17, 25, 18, 26>(data);
    avx2_cas::cas2<T>::template cas<18, 25, 17, 23>(data);
    avx2_cas::cas<T, 18, 24>(data);
    avx2_cas::cas3<T>::template cas<18, 23, 19, 27, 20, 28>(data);
    avx2_cas::cas3<T>::template cas<20, 27, 21, 29, 22, 30>(data);
    avx2_cas::cas2<T>::template cas<22, 29, 21, 27>(data);
    avx2_cas::cas<T, 22, 28>(data);
    avx2_cas::cas3<T>::template cas<22, 27, 19, 23, 20, 24>(data);
    avx2_cas::cas3<T>::template cas<20, 23, 21, 25, 22, 26>(data);
    avx2_cas::cas2<T>::template cas<22, 25, 21, 23>(data);
    avx2_cas::cas<T, 22, 24>(data);
    avx2_cas::cas2<T>::template cas<22, 23, 0, 16>(data);
    avx2_cas::cas3<T>::template cas<0, 15, 1, 17, 2, 18>(data);
    avx2_cas::cas2<T>::template cas<2, 17, 1, 15>(data);
    avx2_cas::cas<T, 2, 16>(data);
    avx2_cas::cas3<T>::template cas<2, 15, 3, 19, 4, 20>(data);
    avx2_cas::cas3<T>::template cas<4, 19, 5, 21, 6, 22>(data);
    avx2_cas::cas2<T>::template cas<6, 21, 5, 19>(data);
    avx2_cas::cas<T, 6, 20>(data);
    avx2_cas::cas3<T>::template cas<6, 19, 3, 15, 4, 16>(data);
    avx2_cas::cas3<T>::template cas<4, 15, 5, 17, 6, 18>(data);
    avx2_cas::cas2<T>::template cas<6, 17, 5, 15>(data);
    avx2_cas::cas<T, 6, 16>(data);
    avx2_cas::cas3<T>::template cas<6, 15, 7, 23, 8, 24>(data);
    avx2_cas::cas3<T>::template cas<8, 23, 9, 25, 10, 26>(data);
    avx2_cas::cas2<T>::template cas<10, 25, 9, 23>(data);
    avx2_cas::cas<T, 10, 24>(data);
    avx2_cas::cas3<T>::template cas<10, 23, 11, 27, 12, 28>(data);
    avx2_cas::cas3<T>::template cas<12, 27, 13, 29, 14, 30>(data);
    avx2_cas::cas2<T>::template cas<14, 29, 13, 27>(data);
    avx2_cas::cas<T, 14, 28>(data);
    avx2_cas::cas3<T>::template cas<14, 27, 11, 23, 12, 24>(data);
    avx2_cas::cas3<T>::template cas<12, 23, 13, 25, 14, 26>(data);
    avx2_cas::cas2<T>::template cas<14, 25, 13, 23>(data);
    avx2_cas::cas<T, 14, 24>(data);
    avx2_cas::cas3<T>::template cas<14, 23, 7, 15, 8, 16>(data);
    avx2_cas::cas3<T>::template cas<8, 15, 9, 17, 10, 18>(data);
    avx2_cas::cas2<T>::template cas<10, 17, 9, 15>(data);
    avx2_cas::cas<T, 10, 16>(data);
    avx2_cas::cas3<T>::template cas<10, 15, 11, 19, 12, 20>(data);
    avx2_cas::cas3<T>::template cas<12, 19, 13, 21, 14, 22>(data);
    avx2_cas::cas2<T>::template cas<14, 21, 13, 19>(data);
    avx2_cas::cas<T, 14, 20>(data);
    avx2_cas::cas3<T>::template cas<14, 19, 11, 15, 12, 16>(data);
    avx2_cas::cas3<T>::template cas<12, 15, 13, 17, 14, 18>(data);
    avx2_cas::cas2<T>::template cas<14, 17, 13, 15>(data);
    avx2_cas::cas<T, 14, 16>(data);
    avx2_cas::cas<T, 14, 15>(data);
}

template <>
template <typename T>
inline
void
vector_sort::internal::net<32>::sort(T * const data)
{
    avx2_cas::cas2<T>::template cas<0, 1, 2, 3>(data);
    avx2_cas::cas2<T>::template cas<0, 2, 1, 3>(data);
    avx2_cas::cas3<T>::template cas<1, 2, 4, 5, 6, 7>(data);
    avx2_cas::cas2<T>::template cas<4, 6, 5, 7>(data);
    avx2_cas::cas2<T>::template cas<5, 6, 0, 4>(data);
    avx2_cas::cas<T, 1, 5>(data);
    avx2_cas::cas3<T>::template cas<1, 4, 2, 6, 3, 7>(data);
    avx2_cas::cas2<T>::template cas<3, 6, 2, 4>(data);
    avx2_cas::cas<T, 3, 5>(data);
    avx2_cas::cas3<T>::template cas<3, 4, 8, 9, 10, 11>(data);
    avx2_cas::cas2<T>::template cas<8, 10, 9, 11>(data);
    avx2_cas::cas3<T>::template cas<9, 10, 12, 13, 14, 15>(data);
    avx2_cas::cas2<T>::template cas<12, 14, 13, 15>(data);
    avx2_cas::cas2<T>::template cas<13, 14, 8, 12>(data);
    avx2_cas::cas<T, 9, 13>(data);
    avx2_cas::cas3<T>::template cas<9, 12, 10, 14, 11, 15>(data);
    avx2_cas::cas2<T>::template cas<11, 14, 10, 12>(data);
    avx2_cas::cas<T, 11, 13>(data);
    avx2_cas::cas3<T>::template cas<11, 12, 0, 8, 1, 9>(data);
    avx2_cas::cas3<T>::template cas<1, 8, 2, 10, 3, 11>(data);
    avx2_cas::cas2<T>::template cas<3, 10, 2, 8>(data);
    avx2_cas::cas<T, 3, 9>(data);
    avx2_cas::cas3<T>::template cas<3, 8, 4, 12, 5, 13>(data);
    avx2_cas::cas3<T>::template cas<5, 12, 6, 14, 7, 15>(data);
    avx2_cas::cas2<T>::template cas<7, 14, 6, 12>(data);
    avx2_cas::cas<T, 7, 13>(data);
    avx2_cas::cas3<T>::template cas<7, 12, 4, 8, 5, 9>(data);
    avx2_cas::cas3<T>::template cas<5, 8, 6, 10, 7, 11>(data);
    avx2_cas::cas2<T>::template cas<7, 10, 6, 8>(data);
    avx2_cas::cas<T, 7, 9>(data);
    avx2_cas::cas3<T>::template cas<7, 8, 16, 17, 18, 19>(data);
    avx2_cas::cas2<T>::template cas<16, 18, 17, 19>(data);
    avx2_cas::cas3<T>::template cas<17, 18, 20, 21, 22, 23>(data);
    avx2_cas::cas2<T>::template cas<20, 22, 21, 23>(data);
    avx2_cas::cas2<T>::template cas<21, 22, 16, 20>(data);
    avx2_cas::cas<T, 17, 21>(data);
    avx2_cas::cas3<T>::template cas<17, 20, 18, 22, 19, 23>(data);
    avx2_cas::cas2<T>::template cas<19, 22, 18, 20>(data);
    avx2_cas::cas<T, 19, 21>(data);
    avx2_cas::cas3<T>::template cas<19, 20, 24, 25, 26, 27>(data);
    avx2_cas::cas2<T>::template cas<24, 26, 25, 27>(data);
    avx2_cas::cas3<T>::template cas<25, 26, 28, 29, 30, 31>(data);
    avx2_cas::cas2<T>::template cas<28, 30, 29, 31>(data);
    avx2_cas::cas2<T>::template cas<29, 30, 24, 28>(data);
    avx2_cas::cas<T, 25, 29>(data);
    avx2_cas::cas3<T>::template cas<25, 28, 26, 30, 27, 31>(data);
    avx2_cas::cas2<T>::template cas<27, 30, 26, 28>(data);
    avx2_cas::cas<T, 27, 29>(data);
    avx2_cas::cas3<T>::template cas<27, 28, 16, 24, 17, 25>(data);
    avx2_cas::cas3<T>::template cas<17, 24, 18, 26, 19, 27>(data);
    avx2_cas::cas2<T>::template cas<19, 26, 18, 24>(data);
    avx2_cas::cas<T, 19, 25>(data);
    avx2_cas::cas3<T>::template cas<19, 24, 20, 28, 21, 29>(data);
    avx2_cas::cas3<T>::template cas<21, 28, 22, 30, 23, 31>(data);
    avx2_cas::cas2<T>::template cas<23, 30, 22, 28>(data);
    avx2_cas::cas<T, 23, 29>(data);
    avx2_cas::cas3<T>::template cas<23, 28, 20, 24, 21, 25>(data);
    avx2_cas::cas3<T>::template cas<21, 24, 22, 26, 23, 27>(data);
    avx2_cas::cas2<T>::template cas<23, 26, 22, 24>(data);
    avx2_cas::cas<T, 23, 25>(data);
    avx2_cas::cas3<T>::template cas<23, 24, 0, 16, 1, 17>(data);
    avx2_cas::cas3<T>::template cas<1, 16, 2, 18, 3, 19>(data);
    avx2_cas::cas2<T>::template cas<3, 18, 2, 16>(data);
    avx2_cas::cas<T, 3, 17>(data);
    avx2_cas::cas3<T>::template cas<3, 16, 4, 20, 5, 21>(data);
    avx2_cas::cas3<T>::template cas<5, 20, 6, 22, 7, 23>(data);
    avx2_cas::cas2<T>::template cas<7, 22, 6, 20>(data);
    avx2_cas::cas<T, 7, 21>(data);
    avx2_cas::cas3<T>::template cas<7, 20, 4, 16, 5, 17>(data);
    avx2_cas::cas3<T>::template cas<5, 16, 6, 18, 7, 19>(data);
    avx2_cas::cas2<T>::template cas<7, 18, 6, 16>(data);
    avx2_cas::cas<T, 7, 17>(data);
    avx2_cas::cas3<T>::template cas<7, 16, 8, 24, 9, 25>(data);
    avx2_cas::cas3<T>::template cas<9, 24, 10, 26, 11, 27>(data);
    avx2_cas::cas2<T>::template cas<11, 26, 10, 24>(data);
    avx2_cas::cas<T, 11, 25>(data);
    avx2_cas::cas3<T>::template cas<11, 24, 12, 28, 13, 29>(data);
    avx2_cas::cas3<T>::template cas<13, 28, 14, 30, 15, 31>(data);
    avx2_cas::cas2<T>::template cas<15, 30, 14, 28>(data);
    avx2_cas::cas<T, 15, 29>(data);
    avx2_cas::cas3<T>::template cas<15, 28, 12, 24, 13, 25>(data);
    avx2_cas::cas3<T>::template cas<13, 24, 14, 26, 15, 27>(data);
    avx2_cas::cas2<T>::template cas<15, 26, 14, 24>(data);
    avx2_cas::cas<T, 15, 25>(data);
    avx2_cas::cas3<T>::template cas<15, 24, 8, 16, 9, 17>(data);
    avx2_cas::cas3<T>::template cas<9, 16, 10, 18, 11, 19>(data);
    avx2_cas::cas2<T>::template cas<11, 18, 10, 16>(data);
    avx2_cas::cas<T, 11, 17>(data);
    avx2_cas::cas3<T>::template cas<11, 16, 12, 20, 13, 21>(data);
    avx2_cas::cas3<T>::template cas<13, 20, 14, 22, 15, 23>(data);
    avx2_cas::cas2<T>::template cas<15, 22, 14, 20>(data);
    avx2_cas::cas<T, 15, 21>(data);
    avx2_cas::cas3<T>::template cas<15, 20, 12, 16, 13, 17>(data);
    avx2_cas::cas3<T>::template cas<13, 16, 14, 18, 15, 19>(data);
    avx2_cas::cas2<T>::template cas<15, 18, 14, 16>(data);
    avx2_cas::cas<T, 15, 17>(data);
    avx2_cas::cas<T, 15, 16>(data);
}

template <typename T>
inline void vector_sort::internal::static_sort(T * const data, size_t const n)
{
    if (n < 2)
    {
        return;
    }
    switch (n)
    {
    case 2:
        vector_sort::internal::net<2>::sort(data);
    break;
    case 3:
        vector_sort::internal::net<3>::sort(data);
    break;
    case 4:
        vector_sort::internal::net<4>::sort(data);
    break;
    case 5:
        vector_sort::internal::net<5>::sort(data);
    break;
    case 6:
        vector_sort::internal::net<6>::sort(data);
    break;
    case 7:
        vector_sort::internal::net<7>::sort(data);
    break;
    case 8:
        vector_sort::internal::net<8>::sort(data);
    break;
    case 9:
        vector_sort::internal::net<9>::sort(data);
    break;
    case 10:
        vector_sort::internal::net<10>::sort(data);
    break;
    case 11:
        vector_sort::internal::net<11>::sort(data);
    break;
    case 12:
        vector_sort::internal::net<12>::sort(data);
    break;
    case 13:
        vector_sort::internal::net<13>::sort(data);
    break;
    case 14:
        vector_sort::internal::net<14>::sort(data);
    break;
    case 15:
        vector_sort::internal::net<15>::sort(data);
    break;
    case 16:
        vector_sort::internal::net<16>::sort(data);
    break;
    case 17:
        vector_sort::internal::net<17>::sort(data);
    break;
    case 18:
        vector_sort::internal::net<18>::sort(data);
    break;
    case 19:
        vector_sort::internal::net<19>::sort(data);
    break;
    case 20:
        vector_sort::internal::net<20>::sort(data);
    break;
    case 21:
        vector_sort::internal::net<21>::sort(data);
    break;
    case 22:
        vector_sort::internal::net<22>::sort(data);
    break;
    case 23:
        vector_sort::internal::net<23>::sort(data);
    break;
    case 24:
        vector_sort::internal::net<24>::sort(data);
    break;
    case 25:
        vector_sort::internal::net<25>::sort(data);
    break;
    case 26:
        vector_sort::internal::net<26>::sort(data);
    break;
    case 27:
        vector_sort::internal::net<27>::sort(data);
    break;
    case 28:
        vector_sort::internal::net<28>::sort(data);
    break;
    case 29:
        vector_sort::internal::net<29>::sort(data);
    break;
    case 30:
        vector_sort::internal::net<30>::sort(data);
    break;
    case 31:
        vector_sort::internal::net<31>::sort(data);
    break;
    case 32:
        vector_sort::internal::net<32>::sort(data);
    break;
    default: assert(false);
    }
}

#endif // SORTING_NETWORK_HPP

/*
      |                |            _)                |                           |   
   _` |  |   |   _` |  |      __ \   | \ \   /  _ \   __|       __|   _ \    __|  __| 
  (   |  |   |  (   |  |      |   |  |  \ \ /  (   |  |       \__ \  (   |  |     |   
 \__,_| \__,_| \__,_| _|      .__/  _|   \_/  \___/  \__|     ____/ \___/  _|    \__| 
                             _|                                                       
                                        Copyright (c) 2018 M. Welsch <michael@welsch.one>
                            Copyright (c) 2018 S. D. Adams <s.d.adams.software@gmail.com>
*/


namespace vector_sort
{
    namespace internal // Implementation details.
    {
        extern bool has_avx2;
        extern bool initialized_avx2;
        void init_avx2();

        static __attribute__((aligned(64))) int PTL[16 * 8] = {
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            2, 3, 0, 0, 0, 0, 0, 0,
            0, 1, 2, 3, 0, 0, 0, 0,
            4, 5, 0, 0, 0, 0, 0, 0,
            0, 1, 4, 5, 0, 0, 0, 0,
            2, 3, 4, 5, 0, 0, 0, 0,
            0, 1, 2, 3, 4, 5, 0, 0,
            6, 7, 0, 0, 0, 0, 0, 0,
            0, 1, 6, 7, 0, 0, 0, 0,
            2, 3, 6, 7, 0, 0, 0, 0,
            0, 1, 2, 3, 6, 7, 0, 0,
            4, 5, 6, 7, 0, 0, 0, 0,
            0, 1, 4, 5, 6, 7, 0, 0,
            2, 3, 4, 5, 6, 7, 0, 0,
            0, 1, 2, 3, 4, 5, 6, 7
        };

        static __attribute__((aligned(64))) int PTL2[256 * 8] = {
            0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,
            1,0,0,0,0,0,0,0,
            0,1,0,0,0,0,0,0,
            2,0,0,0,0,0,0,0,
            0,2,0,0,0,0,0,0,
            1,2,0,0,0,0,0,0,
            0,1,2,0,0,0,0,0,
            3,0,0,0,0,0,0,0,
            0,3,0,0,0,0,0,0,
            1,3,0,0,0,0,0,0,
            0,1,3,0,0,0,0,0,
            2,3,0,0,0,0,0,0,
            0,2,3,0,0,0,0,0,
            1,2,3,0,0,0,0,0,
            0,1,2,3,0,0,0,0,
            4,0,0,0,0,0,0,0,
            0,4,0,0,0,0,0,0,
            1,4,0,0,0,0,0,0,
            0,1,4,0,0,0,0,0,
            2,4,0,0,0,0,0,0,
            0,2,4,0,0,0,0,0,
            1,2,4,0,0,0,0,0,
            0,1,2,4,0,0,0,0,
            3,4,0,0,0,0,0,0,
            0,3,4,0,0,0,0,0,
            1,3,4,0,0,0,0,0,
            0,1,3,4,0,0,0,0,
            2,3,4,0,0,0,0,0,
            0,2,3,4,0,0,0,0,
            1,2,3,4,0,0,0,0,
            0,1,2,3,4,0,0,0,
            5,0,0,0,0,0,0,0,
            0,5,0,0,0,0,0,0,
            1,5,0,0,0,0,0,0,
            0,1,5,0,0,0,0,0,
            2,5,0,0,0,0,0,0,
            0,2,5,0,0,0,0,0,
            1,2,5,0,0,0,0,0,
            0,1,2,5,0,0,0,0,
            3,5,0,0,0,0,0,0,
            0,3,5,0,0,0,0,0,
            1,3,5,0,0,0,0,0,
            0,1,3,5,0,0,0,0,
            2,3,5,0,0,0,0,0,
            0,2,3,5,0,0,0,0,
            1,2,3,5,0,0,0,0,
            0,1,2,3,5,0,0,0,
            4,5,0,0,0,0,0,0,
            0,4,5,0,0,0,0,0,
            1,4,5,0,0,0,0,0,
            0,1,4,5,0,0,0,0,
            2,4,5,0,0,0,0,0,
            0,2,4,5,0,0,0,0,
            1,2,4,5,0,0,0,0,
            0,1,2,4,5,0,0,0,
            3,4,5,0,0,0,0,0,
            0,3,4,5,0,0,0,0,
            1,3,4,5,0,0,0,0,
            0,1,3,4,5,0,0,0,
            2,3,4,5,0,0,0,0,
            0,2,3,4,5,0,0,0,
            1,2,3,4,5,0,0,0,
            0,1,2,3,4,5,0,0,
            6,0,0,0,0,0,0,0,
            0,6,0,0,0,0,0,0,
            1,6,0,0,0,0,0,0,
            0,1,6,0,0,0,0,0,
            2,6,0,0,0,0,0,0,
            0,2,6,0,0,0,0,0,
            1,2,6,0,0,0,0,0,
            0,1,2,6,0,0,0,0,
            3,6,0,0,0,0,0,0,
            0,3,6,0,0,0,0,0,
            1,3,6,0,0,0,0,0,
            0,1,3,6,0,0,0,0,
            2,3,6,0,0,0,0,0,
            0,2,3,6,0,0,0,0,
            1,2,3,6,0,0,0,0,
            0,1,2,3,6,0,0,0,
            4,6,0,0,0,0,0,0,
            0,4,6,0,0,0,0,0,
            1,4,6,0,0,0,0,0,
            0,1,4,6,0,0,0,0,
            2,4,6,0,0,0,0,0,
            0,2,4,6,0,0,0,0,
            1,2,4,6,0,0,0,0,
            0,1,2,4,6,0,0,0,
            3,4,6,0,0,0,0,0,
            0,3,4,6,0,0,0,0,
            1,3,4,6,0,0,0,0,
            0,1,3,4,6,0,0,0,
            2,3,4,6,0,0,0,0,
            0,2,3,4,6,0,0,0,
            1,2,3,4,6,0,0,0,
            0,1,2,3,4,6,0,0,
            5,6,0,0,0,0,0,0,
            0,5,6,0,0,0,0,0,
            1,5,6,0,0,0,0,0,
            0,1,5,6,0,0,0,0,
            2,5,6,0,0,0,0,0,
            0,2,5,6,0,0,0,0,
            1,2,5,6,0,0,0,0,
            0,1,2,5,6,0,0,0,
            3,5,6,0,0,0,0,0,
            0,3,5,6,0,0,0,0,
            1,3,5,6,0,0,0,0,
            0,1,3,5,6,0,0,0,
            2,3,5,6,0,0,0,0,
            0,2,3,5,6,0,0,0,
            1,2,3,5,6,0,0,0,
            0,1,2,3,5,6,0,0,
            4,5,6,0,0,0,0,0,
            0,4,5,6,0,0,0,0,
            1,4,5,6,0,0,0,0,
            0,1,4,5,6,0,0,0,
            2,4,5,6,0,0,0,0,
            0,2,4,5,6,0,0,0,
            1,2,4,5,6,0,0,0,
            0,1,2,4,5,6,0,0,
            3,4,5,6,0,0,0,0,
            0,3,4,5,6,0,0,0,
            1,3,4,5,6,0,0,0,
            0,1,3,4,5,6,0,0,
            2,3,4,5,6,0,0,0,
            0,2,3,4,5,6,0,0,
            1,2,3,4,5,6,0,0,
            0,1,2,3,4,5,6,0,
            7,0,0,0,0,0,0,0,
            0,7,0,0,0,0,0,0,
            1,7,0,0,0,0,0,0,
            0,1,7,0,0,0,0,0,
            2,7,0,0,0,0,0,0,
            0,2,7,0,0,0,0,0,
            1,2,7,0,0,0,0,0,
            0,1,2,7,0,0,0,0,
            3,7,0,0,0,0,0,0,
            0,3,7,0,0,0,0,0,
            1,3,7,0,0,0,0,0,
            0,1,3,7,0,0,0,0,
            2,3,7,0,0,0,0,0,
            0,2,3,7,0,0,0,0,
            1,2,3,7,0,0,0,0,
            0,1,2,3,7,0,0,0,
            4,7,0,0,0,0,0,0,
            0,4,7,0,0,0,0,0,
            1,4,7,0,0,0,0,0,
            0,1,4,7,0,0,0,0,
            2,4,7,0,0,0,0,0,
            0,2,4,7,0,0,0,0,
            1,2,4,7,0,0,0,0,
            0,1,2,4,7,0,0,0,
            3,4,7,0,0,0,0,0,
            0,3,4,7,0,0,0,0,
            1,3,4,7,0,0,0,0,
            0,1,3,4,7,0,0,0,
            2,3,4,7,0,0,0,0,
            0,2,3,4,7,0,0,0,
            1,2,3,4,7,0,0,0,
            0,1,2,3,4,7,0,0,
            5,7,0,0,0,0,0,0,
            0,5,7,0,0,0,0,0,
            1,5,7,0,0,0,0,0,
            0,1,5,7,0,0,0,0,
            2,5,7,0,0,0,0,0,
            0,2,5,7,0,0,0,0,
            1,2,5,7,0,0,0,0,
            0,1,2,5,7,0,0,0,
            3,5,7,0,0,0,0,0,
            0,3,5,7,0,0,0,0,
            1,3,5,7,0,0,0,0,
            0,1,3,5,7,0,0,0,
            2,3,5,7,0,0,0,0,
            0,2,3,5,7,0,0,0,
            1,2,3,5,7,0,0,0,
            0,1,2,3,5,7,0,0,
            4,5,7,0,0,0,0,0,
            0,4,5,7,0,0,0,0,
            1,4,5,7,0,0,0,0,
            0,1,4,5,7,0,0,0,
            2,4,5,7,0,0,0,0,
            0,2,4,5,7,0,0,0,
            1,2,4,5,7,0,0,0,
            0,1,2,4,5,7,0,0,
            3,4,5,7,0,0,0,0,
            0,3,4,5,7,0,0,0,
            1,3,4,5,7,0,0,0,
            0,1,3,4,5,7,0,0,
            2,3,4,5,7,0,0,0,
            0,2,3,4,5,7,0,0,
            1,2,3,4,5,7,0,0,
            0,1,2,3,4,5,7,0,
            6,7,0,0,0,0,0,0,
            0,6,7,0,0,0,0,0,
            1,6,7,0,0,0,0,0,
            0,1,6,7,0,0,0,0,
            2,6,7,0,0,0,0,0,
            0,2,6,7,0,0,0,0,
            1,2,6,7,0,0,0,0,
            0,1,2,6,7,0,0,0,
            3,6,7,0,0,0,0,0,
            0,3,6,7,0,0,0,0,
            1,3,6,7,0,0,0,0,
            0,1,3,6,7,0,0,0,
            2,3,6,7,0,0,0,0,
            0,2,3,6,7,0,0,0,
            1,2,3,6,7,0,0,0,
            0,1,2,3,6,7,0,0,
            4,6,7,0,0,0,0,0,
            0,4,6,7,0,0,0,0,
            1,4,6,7,0,0,0,0,
            0,1,4,6,7,0,0,0,
            2,4,6,7,0,0,0,0,
            0,2,4,6,7,0,0,0,
            1,2,4,6,7,0,0,0,
            0,1,2,4,6,7,0,0,
            3,4,6,7,0,0,0,0,
            0,3,4,6,7,0,0,0,
            1,3,4,6,7,0,0,0,
            0,1,3,4,6,7,0,0,
            2,3,4,6,7,0,0,0,
            0,2,3,4,6,7,0,0,
            1,2,3,4,6,7,0,0,
            0,1,2,3,4,6,7,0,
            5,6,7,0,0,0,0,0,
            0,5,6,7,0,0,0,0,
            1,5,6,7,0,0,0,0,
            0,1,5,6,7,0,0,0,
            2,5,6,7,0,0,0,0,
            0,2,5,6,7,0,0,0,
            1,2,5,6,7,0,0,0,
            0,1,2,5,6,7,0,0,
            3,5,6,7,0,0,0,0,
            0,3,5,6,7,0,0,0,
            1,3,5,6,7,0,0,0,
            0,1,3,5,6,7,0,0,
            2,3,5,6,7,0,0,0,
            0,2,3,5,6,7,0,0,
            1,2,3,5,6,7,0,0,
            0,1,2,3,5,6,7,0,
            4,5,6,7,0,0,0,0,
            0,4,5,6,7,0,0,0,
            1,4,5,6,7,0,0,0,
            0,1,4,5,6,7,0,0,
            2,4,5,6,7,0,0,0,
            0,2,4,5,6,7,0,0,
            1,2,4,5,6,7,0,0,
            0,1,2,4,5,6,7,0,
            3,4,5,6,7,0,0,0,
            0,3,4,5,6,7,0,0,
            1,3,4,5,6,7,0,0,
            0,1,3,4,5,6,7,0,
            2,3,4,5,6,7,0,0,
            0,2,3,4,5,6,7,0,
            1,2,3,4,5,6,7,0,
            0,1,2,3,4,5,6,7
        };

        template <typename T>
        void
        partition(T *in,
                  T *tmp1,
                  T *tmp2,
                  unsigned long n,
                  unsigned long &out_bottom_offset,
                  unsigned long &out_bottom_count,
                  unsigned long &out_middle_offset,
                  unsigned long &out_middle_count,
                  unsigned long &out_top_offset,
                  unsigned long &out_top_count);

        template <typename T>
        void
        insertion_sort(T *d, unsigned long n);

        template <typename T>
        void
        partition_avx2(T *&in,
                       unsigned long &n,
                       T pivot1,
                       T pivot2,
                       T *&bottom,
                       T *&middle,
                       T *&top);

        template <typename T>
        void
        partition_avx2(T *&in,
                       unsigned long &n,
                       T pivot,
                       T *&bottom,
                       T *&middle,
                       T *&top);

        void
        partition_4d(double *&in,
                     __m256d &PIVOT_L,
                     __m256d &PIVOT_H,
                     double *&bottom,
                     double *&middle,
                     double *&top);
        void
        partition_4d(double *&in,
                     __m256d &PIVOT,
                     double *&bottom,
                     double *&middle,
                     double *&top);

        void
        partition_8f(float *&in,
                     __m256 &PIVOT_L,
                     __m256 &PIVOT_H,
                     float *&bottom,
                     float *&middle,
                     float *&top);
        void
        partition_8f(float *&in,
                     __m256 &PIVOT,
                     float *&bottom,
                     float *&middle,
                     float *&top);

        void
        partition_8i(unsigned *&in,
                     __m256i &PIVOT_L,
                     __m256i &PIVOT_H,
                     unsigned *&bottom,
                     unsigned *&middle,
                     unsigned *&top);
        void
        partition_8i(unsigned *&in,
                     __m256i &PIVOT,
                     unsigned *&bottom,
                     unsigned *&middle,
                     unsigned *&top);

        template <typename T>
        void
        sort_recursive(T *unsorted_array, T *tmp_array1, T *tmp_array2, unsigned long n);
    } // end namespace internal
} // end namespace vector_sort

template <typename T>
void vector_sort::internal::insertion_sort(T * const data, unsigned long const count)
{
    if (count <= 1)
    {
        return;
    }

    long i, j;
    for (i = 1; i < count; i++)
    {
        T tmp = data[i];
        for (j = i; j >= 1 && tmp < data[j - 1]; j--)
            data[j] = data[j - 1];
        data[j] = tmp;
    }
}

// Separates 4 doubles into 3 partitions using AVX2.
inline
void
vector_sort::internal::partition_4d(double *&in,
                                    __m256d &PIVOT_L,
                                    __m256d &PIVOT_H,
                                    double *&bottom,
                                    double *&middle,
                                    double *&top)
{
    __m256d INPUT;
    __m256d LOW;
    __m256d MED;
    __m256d HIGH;

    INPUT = _mm256_loadu_pd(in);                     // Load 4 doubles.
    MED = _mm256_cmp_pd(INPUT, PIVOT_L, _CMP_LE_OQ); // Compare less than or equal to the low pivot.
    unsigned long lmask = _mm256_movemask_pd(MED);   // Convert the 4 comparisons to bits in lmask.
    MED = _mm256_cmp_pd(INPUT, PIVOT_H, _CMP_GE_OQ); // Compare greater than or equal to the high pivot.
    unsigned long hmask = _mm256_movemask_pd(MED);   // Convert the 4 comparisons to bits in hmask.
    unsigned long mmask = 0xF & (~lmask & ~hmask);   // Compute the comparisons that produced false for both tests.

    // Partition using the comparisons to index into a permuation array.
    LOW  = (__m256d)_mm256_permutevar8x32_ps((__m256)INPUT, ((__m256i *)PTL)[lmask]);
    MED  = (__m256d)_mm256_permutevar8x32_ps((__m256)INPUT, ((__m256i *)PTL)[mmask]);
    HIGH = (__m256d)_mm256_permutevar8x32_ps((__m256)INPUT, ((__m256i *)PTL)[hmask]);

    // Store unaligned the partitioned values in the corresponding array.
    _mm256_storeu_pd(bottom, LOW);
    _mm256_storeu_pd(middle, MED);
    _mm256_storeu_pd(top, HIGH);

    // Shift the partition bounds by the number of elements stored.
    bottom += __builtin_popcount(lmask);
    middle += __builtin_popcount(mmask);
    top    += __builtin_popcount(hmask);

    in += 4;
}
inline
void
vector_sort::internal::partition_4d(double *&in,
                                    __m256d &PIVOT,
                                    double *&bottom,
                                    double *&middle,
                                    double *&top)
{
    __m256d INPUT;
    __m256d LOW;
    __m256d MED;
    __m256d HIGH;

    INPUT = _mm256_loadu_pd(in);                   // Load 4 doubles.
    MED = _mm256_cmp_pd(INPUT, PIVOT, _CMP_LT_OQ); // Compare less than or equal to the low pivot.
    unsigned long lmask = _mm256_movemask_pd(MED); // Convert the 4 comparisons to bits in lmask.
    MED = _mm256_cmp_pd(INPUT, PIVOT, _CMP_GT_OQ); // Compare greater than or equal to the high pivot.
    unsigned long hmask = _mm256_movemask_pd(MED); // Convert the 4 comparisons to bits in hmask.
    unsigned long mmask = 0xF & (~lmask & ~hmask); // Compute the comparisons that produced false for both tests.

    // Partition using the comparisons to index into a permuation array.
    LOW  = (__m256d)_mm256_permutevar8x32_ps((__m256)INPUT, ((__m256i *)PTL)[lmask]);
    MED  = (__m256d)_mm256_permutevar8x32_ps((__m256)INPUT, ((__m256i *)PTL)[mmask]);
    HIGH = (__m256d)_mm256_permutevar8x32_ps((__m256)INPUT, ((__m256i *)PTL)[hmask]);

    // Store unaligned the partitioned values in the corresponding array.
    _mm256_storeu_pd(bottom, LOW);
    _mm256_storeu_pd(middle, MED);
    _mm256_storeu_pd(top, HIGH);

    // Shift the partition bounds by the number of elements stored.
    bottom += __builtin_popcount(lmask);
    middle += __builtin_popcount(mmask);
    top    += __builtin_popcount(hmask);

    in += 4;
}

// Separates 8 floats into 3 partitions using AVX2.
inline
void
vector_sort::internal::partition_8f(float *&in,
                                    __m256 &PIVOT_L,
                                    __m256 &PIVOT_H,
                                    float *&bottom,
                                    float *&middle,
                                    float *&top)
{
    __m256 INPUT;
    __m256 LOW;
    __m256 MED;
    __m256 HIGH;

    INPUT = _mm256_loadu_ps(in);                    // Load 8 floats.
    MED = _mm256_cmp_ps(INPUT, PIVOT_L, _CMP_LE_OQ);         // Compare less than or equal to the low pivot.
    unsigned long lmask = _mm256_movemask_ps(MED);  // Convert the 8 comparisons to bits in lmask.
    MED = _mm256_cmp_ps(INPUT, PIVOT_H, _CMP_GE_OQ);        // Compare greater than or equal to the high pivot.
    unsigned long hmask = _mm256_movemask_ps(MED);  // Convert the 8 comparisons to bits in hmask.
    unsigned long mmask = 0xFF & (~lmask & ~hmask); // Compute the comparisons that produced false for both tests.

    // Partition using the comparisons to index into a permuation array.
    LOW  = _mm256_permutevar8x32_ps(INPUT, ((__m256i *)PTL2)[lmask]);
    MED  = _mm256_permutevar8x32_ps(INPUT, ((__m256i *)PTL2)[mmask]);
    HIGH = _mm256_permutevar8x32_ps(INPUT, ((__m256i *)PTL2)[hmask]);

    // Store unaligned the partitioned values in the corresponding array.
    _mm256_storeu_ps(bottom, LOW);
    _mm256_storeu_ps(middle, MED);
    _mm256_storeu_ps(top, HIGH);

    // Shift the partition bounds by the number of elements stored.
    bottom += __builtin_popcount(lmask);
    middle += __builtin_popcount(mmask);
    top    += __builtin_popcount(hmask);

    in += 8;
}
inline
void
vector_sort::internal::partition_8f(float *&in,
                                    __m256 &PIVOT,
                                    float *&bottom,
                                    float *&middle,
                                    float *&top)
{
    __m256 INPUT;
    __m256 LOW;
    __m256 MED;
    __m256 HIGH;

    INPUT = _mm256_loadu_ps(in);                    // Load 8 floats.
    MED = _mm256_cmp_ps(INPUT, PIVOT, _CMP_LT_OQ);         // Compare less than or equal to the low pivot.
    unsigned long lmask = _mm256_movemask_ps(MED);  // Convert the 8 comparisons to bits in lmask.
    MED = _mm256_cmp_ps(INPUT, PIVOT, _CMP_GT_OQ);        // Compare greater than or equal to the high pivot.
    unsigned long hmask = _mm256_movemask_ps(MED);  // Convert the 8 comparisons to bits in hmask.
    unsigned long mmask = 0xFF & (~lmask & ~hmask); // Compute the comparisons that produced false for both tests.

    // Partition using the comparisons to index into a permuation array.
    LOW  = _mm256_permutevar8x32_ps(INPUT, ((__m256i *)PTL2)[lmask]);
    MED  = _mm256_permutevar8x32_ps(INPUT, ((__m256i *)PTL2)[mmask]);
    HIGH = _mm256_permutevar8x32_ps(INPUT, ((__m256i *)PTL2)[hmask]);

    // Store unaligned the partitioned values in the corresponding array.
    _mm256_storeu_ps(bottom, LOW);
    _mm256_storeu_ps(middle, MED);
    _mm256_storeu_ps(top, HIGH);

    // Shift the partition bounds by the number of elements stored.
    bottom += __builtin_popcount(lmask);
    middle += __builtin_popcount(mmask);
    top    += __builtin_popcount(hmask);

    in += 8;
}

// Separates 8 ints into 3 partitions using AVX2.
inline
void
vector_sort::internal::partition_8i(unsigned *&in,
                                    __m256i &PIVOT_L,
                                    __m256i &PIVOT_H,
                                    unsigned *&bottom,
                                    unsigned *&middle,
                                    unsigned *&top)
{
    __m256i INPUT;
    __m256i LOW;
    __m256i MED;
    __m256i HIGH;

    INPUT = _mm256_loadu_si256((__m256i *)in);       // Load 8 ints.
    __m256i LMASK_LT = _mm256_cmpgt_epi32(PIVOT_L, INPUT);        // Compare less than the low pivot.
    __m256i LMASK_EQ = _mm256_cmpeq_epi32(PIVOT_L, INPUT);        // Compare equal to the low pivot.
    __m256i HMASK_GT = _mm256_cmpgt_epi32(INPUT, PIVOT_H);        // Compare greater than the high pivot.
    __m256i HMASK_EQ = _mm256_cmpeq_epi32(INPUT, PIVOT_H);        // Compare equal to the high pivot.
    __m256i LMASK = _mm256_or_si256(LMASK_LT, LMASK_EQ);
    __m256i HMASK = _mm256_or_si256(HMASK_GT, HMASK_EQ);
    __m256i TMASK = _mm256_set1_epi32(-1);
    __m256i MMASK = _mm256_andnot_si256(HMASK, _mm256_andnot_si256(LMASK, TMASK));

//    unsigned long int lmask = _mm256_movemask_ps(_mm256_castsi256_ps(LMASK));
//    unsigned long int mmask = _mm256_movemask_ps(_mm256_castsi256_ps(MMASK));
//    unsigned long int hmask = _mm256_movemask_ps(_mm256_castsi256_ps(HMASK));

    unsigned long int lmask = _mm256_movemask_ps((__m256)LMASK);
    unsigned long int mmask = _mm256_movemask_ps((__m256)MMASK);
    unsigned long int hmask = _mm256_movemask_ps((__m256)HMASK);

    // Partition using the comparisons to index into a permuation array.
    LOW  = _mm256_permutevar8x32_epi32(INPUT, ((__m256i *)PTL2)[lmask]);
    MED  = _mm256_permutevar8x32_epi32(INPUT, ((__m256i *)PTL2)[mmask]);
    HIGH = _mm256_permutevar8x32_epi32(INPUT, ((__m256i *)PTL2)[hmask]);

    // Store unaligned the partitioned values in the corresponding array.
    _mm256_storeu_si256((__m256i *)bottom, LOW);
    _mm256_storeu_si256((__m256i *)middle, MED);
    _mm256_storeu_si256((__m256i *)top, HIGH);

    // Shift the partition bounds by the number of elements stored.
    unsigned long int lp = __builtin_popcount(lmask);
    unsigned long int mp = __builtin_popcount(mmask);
    unsigned long int hp = __builtin_popcount(hmask);
    assert (lp + mp + hp == 8);
    bottom += lp;
    middle += mp;
    top    += hp;

    in += 8;
}
inline
void
vector_sort::internal::partition_8i(unsigned *&in,
                                    __m256i &PIVOT,
                                    unsigned *&bottom,
                                    unsigned *&middle,
                                    unsigned *&top)
{
    __m256i INPUT;
    __m256i LOW;
    __m256i MED;
    __m256i HIGH;

    INPUT = _mm256_loadu_si256((__m256i *)in);       // Load 8 ints.
    __m256i LMASK = _mm256_cmpgt_epi32(PIVOT, INPUT);        // Compare less than the low pivot.
    __m256i HMASK = _mm256_cmpgt_epi32(INPUT, PIVOT);        // Compare greater than the high pivot.
    __m256i TMASK = _mm256_set1_epi32(-1);
    __m256i MMASK = _mm256_andnot_si256(HMASK, _mm256_andnot_si256(LMASK, TMASK));

//    unsigned long int lmask = _mm256_movemask_ps(_mm256_castsi256_ps(LMASK));
//    unsigned long int mmask = _mm256_movemask_ps(_mm256_castsi256_ps(MMASK));
//    unsigned long int hmask = _mm256_movemask_ps(_mm256_castsi256_ps(HMASK));

    unsigned long int lmask = _mm256_movemask_ps((__m256)LMASK);
    unsigned long int mmask = _mm256_movemask_ps((__m256)MMASK);
    unsigned long int hmask = _mm256_movemask_ps((__m256)HMASK);

    // Partition using the comparisons to index into a permuation array.
    LOW  = _mm256_permutevar8x32_epi32(INPUT, ((__m256i *)PTL2)[lmask]);
    MED  = _mm256_permutevar8x32_epi32(INPUT, ((__m256i *)PTL2)[mmask]);
    HIGH = _mm256_permutevar8x32_epi32(INPUT, ((__m256i *)PTL2)[hmask]);

    // Store unaligned the partitioned values in the corresponding array.
    _mm256_storeu_si256((__m256i *)bottom, LOW);
    _mm256_storeu_si256((__m256i *)middle, MED);
    _mm256_storeu_si256((__m256i *)top, HIGH);

    // Shift the partition bounds by the number of elements stored.
    unsigned long int lp = __builtin_popcount(lmask);
    unsigned long int mp = __builtin_popcount(mmask);
    unsigned long int hp = __builtin_popcount(hmask);
    assert (lp + mp + hp == 8);
    bottom += lp;
    middle += mp;
    top    += hp;

    in += 8;
}

template <>
inline
void
vector_sort::internal::partition_avx2<double>(double *&in,
                                              unsigned long &n,
                                              double pivot1,
                                              double pivot2,
                                              double *&bottom,
                                              double *&middle,
                                              double *&top)
{
    // Duplicate the pivot values 4 times into vector equivalents.
    __m256d PIVOT1 = _mm256_broadcast_sd(&pivot1);
    __m256d PIVOT2 = _mm256_broadcast_sd(&pivot2);

    // Partition the input one 4-element AVX2 vector per loop iteration.
    while (n >= 4)
    {
        partition_4d(in, PIVOT1, PIVOT2, bottom, middle, top);
        n -= 4;
    }
}

template <>
inline
void
vector_sort::internal::partition_avx2<double>(double *&in,
                                              unsigned long &n,
                                              double pivot,
                                              double *&bottom,
                                              double *&middle,
                                              double *&top)
{
    // Duplicate the pivot values 4 times into a vector equivalent.
    __m256d PIVOT = _mm256_broadcast_sd(&pivot);

    // Partition the input one 4-element AVX2 vector per loop iteration.
    while (n >= 4)
    {
        partition_4d(in, PIVOT, bottom, middle, top);
        n -= 4;
    }
}

template <>
inline
void
vector_sort::internal::partition_avx2<float>(float *&in,
                                             unsigned long &n,
                                             float pivot1,
                                             float pivot2,
                                             float *&bottom,
                                             float *&middle,
                                             float *&top)
{
    // Duplicate the pivot values 8 times into vector equivalents.
    __m256 PIVOT1 = _mm256_broadcast_ss(&pivot1);
    __m256 PIVOT2 = _mm256_broadcast_ss(&pivot2);

    // Partition the input one 8-element AVX2 vector per loop iteration.
    while (n >= 8)
    {
        partition_8f(in, PIVOT1, PIVOT2, bottom, middle, top);
        n -= 8;
    }
}

template <>
inline
void
vector_sort::internal::partition_avx2<float>(float *&in,
                                             unsigned long &n,
                                             float pivot,
                                             float *&bottom,
                                             float *&middle,
                                             float *&top)
{
    // Duplicate the pivot values 8 times into vector equivalents.
    __m256 PIVOT = _mm256_broadcast_ss(&pivot);

    // Partition the input one 8-element AVX2 vector per loop iteration.
    while (n >= 8)
    {
        partition_8f(in, PIVOT, bottom, middle, top);
        n -= 8;
    }
}

template <>
inline
void
vector_sort::internal::partition_avx2<unsigned>(unsigned *&in,
                                           unsigned long &n,
                                           unsigned pivot1,
                                           unsigned pivot2,
                                           unsigned *&bottom,
                                           unsigned *&middle,
                                           unsigned *&top)
{
    // Duplicate the pivot values 8 times into vector equivalents.
    __m256i PIVOT1 = _mm256_set1_epi32(pivot1);
    __m256i PIVOT2 = _mm256_set1_epi32(pivot2);

    // Partition the input one 8-element AVX2 vector per loop iteration.
    while (n >= 8)
    {
        partition_8i(in, PIVOT1, PIVOT2, bottom, middle, top);
        n -= 8;
    }
}

template <>
inline
void
vector_sort::internal::partition_avx2<unsigned>(unsigned *&in,
                                           unsigned long &n,
                                           unsigned pivot,
                                           unsigned *&bottom,
                                           unsigned *&middle,
                                           unsigned *&top)
{
    // Duplicate the pivot values 8 times into vector equivalents.
    __m256i PIVOT = _mm256_set1_epi32(pivot);

    // Partition the input one 8-element AVX2 vector per loop iteration.
    while (n >= 8)
    {
        partition_8i(in, PIVOT, bottom, middle, top);
        n -= 8;
    }
}

template <typename T>
void
vector_sort::internal::partition(T *in,
                                 T *tmp1,
                                 T *tmp2,
                                 unsigned long n,
                                 unsigned long &out_bottom_offset,
                                 unsigned long &out_bottom_count,
                                 unsigned long &out_middle_offset,
                                 unsigned long &out_middle_count,
                                 unsigned long &out_top_offset,
                                 unsigned long &out_top_count)
{
    T * const orig_in = in;

    T *bottom = in;
    T *middle = tmp1;
    T *top = tmp2;

    // Compute the two pivots by looking at 5 elements.
    size_t
        quarter = n / 4,
        i1 = 0,
        i5 = n - 1,
        i3 = n / 2,
        i2 = i3 - quarter,
        i4 = i3 + quarter;
    T
        &e1 = in[i1],
        &e2 = in[i2],
        &e3 = in[i3],
        &e4 = in[i4],
        &e5 = in[i5];

    T t; // A temporary to assist swapping.

    // Sort the five elements.
    if (e1 > e2)
        t = e1, e1 = e2, e2 = t;
    if (e4 > e5)
        t = e4, e4 = e5, e5 = t;
    if (e1 > e3)
        t = e1, e1 = e3, e3 = t;
    if (e2 > e3)
        t = e2, e2 = e3, e3 = t;
    if (e1 > e4)
        t = e1, e1 = e4, e4 = t;
    if (e3 > e4)
        t = e3, e3 = e4, e4 = t;
    if (e2 > e5)
        t = e2, e2 = e5, e5 = t;
    if (e2 > e3)
        t = e2, e2 = e3, e3 = t;
    if (e4 > e5)
        t = e4, e4 = e5, e5 = t;

    // Establish the pivot values.
    T pivot1 = e2, pivot2 = e4;

    bool const single_pivoted = (pivot1 == pivot2);

    // The bulk of the partitioning is passed to a type-specific implementation accelerated using AVX2.
    if (!single_pivoted)
    {
        partition_avx2<T>(in, n, pivot1, pivot2, bottom, middle, top);
    }
    else
    {
        partition_avx2<T>(in, n, pivot1, bottom, middle, top);
    }

    // Handle the remainder of elements that can't be vectorized.
    if (!single_pivoted)
    {
        while (n > 0)
        {
            if (in[0] <= pivot1)
            {
                bottom[0] = in[0];
                bottom += 1;
                in += 1;
                n--;
            }
            else if (in[0] >= pivot2)
            {
                top[0] = in[0];
                top += 1;
                in += 1;
                n--;
            }
            else
            {
                middle[0] = in[0];
                middle += 1;
                in += 1;
                n--;
            }
        }
    }
    else
    {
        while (n > 0)
        {
            if (in[0] < pivot1)
            {
                bottom[0] = in[0];
                bottom += 1;
                in += 1;
                n--;
            }
            else if (in[0] > pivot2)
            {
                top[0] = in[0];
                top += 1;
                in += 1;
                n--;
            }
            else
            {
                middle[0] = in[0];
                middle += 1;
                in += 1;
                n--;
            }
        }
    }

    // Rejoin the partitioned elements to the input array.
    unsigned long const bottom_count = ((uint64_t)bottom - (uint64_t)orig_in) / sizeof(T);
    unsigned long const middle_count = ((uint64_t)middle - (uint64_t)tmp1) / sizeof(T);
    unsigned long const top_count    = ((uint64_t)top - (uint64_t)tmp2) / sizeof(T);
    if (middle_count > 0)
    {
        memcpy(bottom, tmp1, middle_count * sizeof(T));
    }
    if (top_count > 0)
    {
        memcpy(bottom + middle_count, tmp2, top_count * sizeof(T));
    }

    out_bottom_offset = 0;
    out_bottom_count = bottom_count;
    out_middle_offset = bottom_count;
    out_middle_count = single_pivoted ? 0 : middle_count;
    out_top_offset = bottom_count + middle_count;
    out_top_count = top_count;
    while (out_bottom_count > 0 && orig_in[out_bottom_count - 1] == pivot1)
    {
        --out_bottom_count;
    }
    while (out_top_count > 0 && orig_in[out_top_offset] == pivot2)
    {
        ++out_top_offset;
        --out_top_count;
    }
    while (out_middle_count > 0 && orig_in[out_middle_offset] == pivot1)
    {
        ++out_middle_offset;
        --out_middle_count;
    }
    while (out_middle_count > 0 && orig_in[out_middle_count - 1] == pivot2)
    {
        --out_middle_count;
    }
}








namespace vector_sort
{
    template <typename T>
    void
    sort(T *data, size_t const count); // Entry-point to this module for plain arrays.

    template <typename T>
    void
    sort(std::vector<T> &vec);         // Entry-point to this module for std::vectors.

    bool has_avx2();


#define VECTOR_SORT_FALLBACK_INSERTION_SORT 0

} // namepace vector_sort

template <typename T>
void
vector_sort::internal::sort_recursive(T *unsorted_array, T *tmp_array1, T *tmp_array2, unsigned long n)
{
    if (n < 2)
    {
        return;
    }

    // Recursively sort until partitions are switchPoint-or-less elements in size.
    unsigned long
            low_offset,
            low_count,
            middle_offset,
            middle_count,
            high_offset,
            high_count;
    while (n > vector_sort::internal::switchPoint)
    {
        vector_sort::internal::partition(unsorted_array,
                                                                tmp_array1,
                                                                tmp_array2,
                                                                n,
                                                                low_offset,
                                                                low_count,
                                                                middle_offset,
                                                                middle_count,
                                                                high_offset,
                                                                high_count);
        n = low_count;
        vector_sort::internal::sort_recursive(&unsorted_array[middle_offset], tmp_array1, tmp_array2, middle_count);
        vector_sort::internal::sort_recursive(&unsorted_array[high_offset],   tmp_array1, tmp_array2, high_count);
    }

#if VECTOR_SORT_FALLBACK_INSERTION_SORT
    // Delegate the sorting of partitions of at most switchPoint elements to insertion sort.
    vector_sort::internal::insertion_sort(unsorted_array, n);
#else
    vector_sort::internal::static_sort<T>(unsorted_array, n);
#endif
}

template <typename T>
void
vector_sort::sort(T * const data, size_t const count)
{
    if (vector_sort::has_avx2())
    {
        T * const tmp1 = (T *)malloc(sizeof(T) * count);
        T * const tmp2 = (T *)malloc(sizeof(T) * count);
    
        vector_sort::internal::sort_recursive(data, tmp1, tmp2, count);
    
        free(tmp1);
        free(tmp2);
    }
    else
    {
        std::sort(data, data + count);
    }
}

template <typename T>
void
vector_sort::sort(std::vector<T> &vec)
{
    T *data = &vec[0];
    size_t const count = vec.size();

    vector_sort::sort(data, count);
}

bool vector_sort::has_avx2()
{
    if (!vector_sort::internal::initialized_avx2)
    {
        vector_sort::internal::init_avx2();
    }
    return vector_sort::internal::has_avx2;
}

bool vector_sort::internal::has_avx2(false);

bool vector_sort::internal::initialized_avx2(false);

void vector_sort::internal::init_avx2()
{
    if (vector_sort::internal::initialized_avx2)
    {
        return;
    }
#if defined(__INTEL_COMPILER)
    if (_may_i_use_cpu_feature(_FEATURE_AVX2))
    {
        vector_sort::internal::has_avx2 = true;
    }
#else
#if defined(__GNUC__) && !defined(__clang__)
    __builtin_cpu_init();
#endif
    if (__builtin_cpu_supports("avx2"))
    {
        vector_sort::internal::has_avx2 = true;
    }
#endif
    vector_sort::internal::initialized_avx2 = true;
}

#endif //SIMD_SORT_H_GUARD