

#ifndef SCORE_VECTOR_H_
#define SCORE_VECTOR_H_

#ifdef __SSSE3__
#include <tmmintrin.h>
#endif

template<typename _score>
struct score_traits
{
	static const unsigned channels = 1;
	enum { zero = 0, byte_size = 4 };
	typedef bool Mask;
};

template<>
struct score_traits<uint8_t>
{
	//static const unsigned channels = 16;
	enum { channels = 16, zero = 0x00, byte_size = 1 };
	typedef uint16_t Mask;
};

template<typename _score>
struct score_vector
{ };

template<>
class score_vector<uint8_t>
{
	public:
	score_vector()
	{
		data_ = _mm_set1_epi8(score_traits<uint8_t>::zero);
	}

	explicit score_vector(int x):
		data_ (_mm_set(x))
	{ }

	explicit score_vector(__m128i data):
		data_ (data)
	{ }

	explicit score_vector(unsigned a, const __m128i &seq)
	{
// 		if(VATParameters::have_ssse3) {
// #ifdef __SSSE3__
// 			set_ssse3(a, seq);
// #else
// 			set_generic(a, seq);
// #endif
// 		} else
			set_generic(a, seq);
	}

// 	void set_ssse3(unsigned a, const __m128i &seq)
// 	{
// #ifdef __SSSE3__
// 		const __m128i *row = reinterpret_cast<const __m128i*>(&score_matrix::get().matrix8u()[a << 5]);

// 		__m128i high_mask = _mm_slli_epi16(_mm_and_si128(seq, _mm_set1_epi8(0x10)), 3);
// 		__m128i seq_low = _mm_or_si128(seq, high_mask);
// 		__m128i seq_high = _mm_or_si128(seq, _mm_xor_si128(high_mask, _mm_set1_epi8(0x80)));

// 		__m128i r1 = _mm_load_si128(row);
// 		__m128i r2 = _mm_load_si128(row+1);
// 		__m128i s1 = _mm_shuffle_epi8(r1, seq_low);
// 		__m128i s2 = _mm_shuffle_epi8(r2, seq_high);
// 		data_ = _mm_or_si128(s1, s2);
// #endif
// 	}

	void set_generic(unsigned a, const __m128i &seq)
	{
		const uint8_t* row (&ScoreMatrix::get().matrix8u()[a<<5]);
		const uint8_t* seq_ptr (reinterpret_cast<const uint8_t*>(&seq));
		uint8_t* dest (reinterpret_cast<uint8_t*>(&data_));
		for(unsigned i=0;i<16;i++)
			*(dest++) = row[*(seq_ptr++)];
	}

	score_vector(const uint8_t* s):
		data_ (_mm_loadu_si128(reinterpret_cast<const __m128i*>(s)))
	{ }

	score_vector operator+(const score_vector &rhs) const
	{
		return score_vector (_mm_adds_epu8(data_, rhs.data_));
	}

	score_vector operator-(const score_vector &rhs) const
	{
		return score_vector (_mm_subs_epu8(data_, rhs.data_));
	}

	score_vector& operator-=(const score_vector &rhs)
	{
		data_ = _mm_subs_epu8(data_, rhs.data_);
		return *this;
	}

	void unbias(const score_vector &bias)
	{ this->operator -=(bias); }

	int operator [](unsigned i) const
	{
		return *(((uint8_t*)&data_)+i);
	}

	void set(unsigned i, uint8_t v)
	{
		*(((uint8_t*)&data_)+i) = v;
	}

	score_vector& max(const score_vector &rhs)
	{
		data_ = _mm_max_epu8(data_, rhs.data_);
		return *this;
	}

	score_vector& min(const score_vector &rhs)
	{
		data_ = _mm_min_epu8(data_, rhs.data_);
		return *this;
	}

	friend score_vector max(const score_vector& lhs, const score_vector &rhs)
	{
		return score_vector (_mm_max_epu8(lhs.data_, rhs.data_));
	}

	friend score_vector min(const score_vector& lhs, const score_vector &rhs)
	{
		return score_vector (_mm_min_epu8(lhs.data_, rhs.data_));
	}

	uint16_t cmpeq(const score_vector &rhs) const
	{
		return _mm_movemask_epi8(_mm_cmpeq_epi8(data_, rhs.data_));
	}

	__m128i cmpeq2(const score_vector &rhs) const
	{
		return _mm_cmpeq_epi8(data_, rhs.data_);
	}

	uint16_t cmpgt(const score_vector &rhs) const
	{
		return _mm_movemask_epi8(_mm_cmpgt_epi8(data_, rhs.data_));
	}

	__m128i data_;

};

#endif /* SCORE_VECTOR_H_ */
