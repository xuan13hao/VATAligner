

#ifndef PACKED_SEQUENCE_H_
#define PACKED_SEQUENCE_H_

#include <vector>
#include "../tools/binary_buffer.h"

using std::vector;

bool has_n(const sequence<const DNA> &seq)
{
	for(unsigned i=0;i<seq.length();++i)
		if(seq[i] == AlphabetAttributes<DNA>::MASK_CHAR)
			return true;
	return false;
}

struct Packed_sequence
{

	// Packed_sequence(const sequence<const Nucleotide> &seq):
	// 	has_n_ (::has_n(seq))
	// {
	// 	if(has_n_)
	// 		pack<Nucleotide,3>(seq);
	// 	else
	// 		pack<Nucleotide,2>(seq);
	// }
	Packed_sequence(const sequence<const DNA> &seq):
		has_n_ (false)
	{ pack<DNA,5>(seq); }
	Packed_sequence(const sequence<const Protein> &seq):
		has_n_ (false)
	{ pack<Protein,5>(seq); }

	Packed_sequence(Binary_buffer::Iterator &it, unsigned len, bool has_n, unsigned b):
		has_n_ (has_n)
	{
		const size_t l = (len*b+7)/8;
		it.read(data_, l);
	}

	template<typename _val>
	void unpack(vector<_val> &dst, unsigned b, unsigned len)
	{
		dst.clear();
		unsigned x = 0, n = 0, l = 0;
		const unsigned mask = (1<<b)-1;
		for(unsigned i=0;i<data_.size();++i) {
			x |= (unsigned)data_[i] << n;
			n += 8;
			while(n >= b && l < len) {
				dst.push_back(x & mask);
				n -= b;
				x >>= b;
				++l;
			}
		}
	}

	const vector<uint8_t>& data() const
	{ return data_; }

	bool has_n() const
	{ return has_n_; }

private:

	template<typename _val, unsigned _b>
	void pack(const sequence<const _val> &seq)
	{
		unsigned x = 0, n = 0;
		for(unsigned i=0;i<seq.length();++i) {
			x |= (unsigned)seq[i] << n;
			n += _b;
			if(n >= 8) {
				data_.push_back(x & 0xff);
				n -= 8;
				x >>= 8;
			}
		}
		if(n > 0)
			data_.push_back(x & 0xff);
	}

	bool has_n_;
	vector<uint8_t> data_;

};

#endif /* PACKED_SEQUENCE_H_ */
