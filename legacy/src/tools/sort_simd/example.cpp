#include <iostream>
#include <emmintrin.h> // for SSE2 intrinsics

int main() {
    float a[4] = { 1.0f, 2.0f, 3.0f, 4.0f };
    float b[4] = { 5.0f, 6.0f, 7.0f, 8.0f };
    float c[4];

    // Load the floats from arrays a and b into two SSE2 registers
    __m128 reg1 = _mm_loadu_ps(a);
    __m128 reg2 = _mm_loadu_ps(b);

    // Add the floats in the two SSE2 registers
    __m128 reg3 = _mm_add_ps(reg1, reg2);

    // Store the result back into array c
    _mm_storeu_ps(c, reg3);

    // Print the result
    std::cout << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << std::endl;

    return 0;
}