
/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef SCORE_PROFILE_H_
#define SCORE_PROFILE_H_

#include <vector>
#include "../Commons/sequence.h"
#include "score_vector.h"

using std::vector;

struct sequence_stream
{
	sequence_stream():
		next (buffer_len),
		mask (0)
	{ }
	inline void reset()
	{
		next = buffer_len;
		mask = 0;
	}
	template<typename _val, typename _score>
	inline const __m128i& get(const typename vector<sequence<const _val> >::const_iterator &begin,
					   const typename vector<sequence<const _val> >::const_iterator &end,
					   unsigned pos,
					   const _score&)
	{
		if(next == buffer_len)
			fill<_val,_score>(begin, end, pos);
		return data_[next++];
	}
	template<typename _val, typename _score>
	inline void fill(const typename vector<sequence<const _val> >::const_iterator &begin,
		  	  const typename vector<sequence<const _val> >::const_iterator &end,
		 	  unsigned pos)
	{
		memset(data_, AlphabetAttributes<_val>::MASK_CHAR, buffer_len*16);
		unsigned n = 0;
		typename vector<sequence<const _val> >::const_iterator it (begin);
		assert(pos < it->length());
		const unsigned read_len (std::min(unsigned(buffer_len), static_cast<unsigned>(it->length())-pos));
		while(it < end) {
			const uint8_t *src (reinterpret_cast<const uint8_t*>(it->data()) + pos);
			_score *dest (reinterpret_cast<_score*>(data_) + n);
			int clip (int(pos) - it->clipping_offset_);
			if((mask & (1 << n)) == 0) {
				if(copy_char(src, dest, mask, n, clip))
				if(read_len > 1 && copy_char(src, dest, mask, n, clip))
				if(read_len > 2 && copy_char(src, dest, mask, n, clip))
				if(read_len > 3) copy_char(src, dest, mask, n, clip);
			}
			++it;
			++n;
		}
		next = 0;
	}
	template<typename _score>
	static inline bool copy_char(const uint8_t*& src, _score*& dest, unsigned &mask, unsigned n, int &clip)
	{
		if(clip++ < 0) {
			dest += 16/sizeof(_score);
			++src;
			return true;
		}
		if(*src == 0xff) {
			mask |= 1 << n;
			return false;
		}
		*dest = *(src++) & 0x7f;
		dest += 16/sizeof(_score);
		return true;
	}
	static const unsigned buffer_len = 4;
	__m128i data_[buffer_len];
	unsigned next;
	unsigned mask;
};

template<typename _score>
struct score_profile
{

	template<typename _val>
	inline void set_(const __m128i &seq)
	{
		assert(sizeof(data_)/sizeof(score_vector<_score>) >= AlphabetAttributes<_val>::ALPHABET_SIZE);
		/*unsigned j = 0;
		do {
			data_[j] = score_vector<_score> (j, seq);
			++j;
			data_[j] = score_vector<_score> (j, seq);
			++j;
			data_[j] = score_vector<_score> (j, seq);
			++j;
			data_[j] = score_vector<_score> (j, seq);
			++j;
		} while(j<24);
		data_[j] = score_vector<_score> (j, seq);
		assert(j+1 == Value_traits<_val>::ALPHABET_SIZE);*/
		for(unsigned j=0;j<AlphabetAttributes<_val>::ALPHABET_SIZE;++j)
			data_[j] = score_vector<_score> (j, seq);
	}

	template<typename _val>
	inline const score_vector<_score>& get(_val i) const
	{ return data_[(int)i]; }

	score_vector<_score> data_[25];

};

#endif /* SCORE_PROFILE_H_ */
