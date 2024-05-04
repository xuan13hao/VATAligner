

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include "../basic/value.h"
#include "../tools/binary_buffer.h"
#include "../tools/text_buffer.h"

template<typename _val>
class sequence
{
	public:
	sequence():
		len_ (0),
		clipping_offset_ (0),
		data_ (0)
	{ }
	sequence(_val *data, size_t len, int clipping_offset = 0):
		len_ (len),
		clipping_offset_ (clipping_offset),
		data_ (data)
	{ }
	size_t length() const
	{
		return len_;
	}
	size_t clipped_length() const
	{
		return len_ - clipping_offset_;
	}
	size_t aligned_clip(unsigned padding) const
	{
		return clipping_offset_ > padding ? clipping_offset_ - padding : 0;
	}
	const _val* data() const
	{
		return data_;
	}
	const _val* clipped_data() const
	{
		return data_ + clipping_offset_;
	}
	const _val* aligned_data(unsigned padding) const
	{
		return data_ + padding;
	}
	const _val& operator [](size_t i) const
	{
		return data_[i];
	}
	_val& operator [](size_t i)
	{
		return data_[i];
	}
	bool empty() const
	{ return len_ == 0; }
	const char* c_str() const
	{ return reinterpret_cast<const char*>(data_); }
	size_t print(char *ptr, unsigned begin, unsigned len) const
	{
		for(unsigned i=begin;i<begin+len;++i)
			*(ptr++) = to_char(data_[i]);
		return len;
	}
	friend std::ostream& operator<<(std::ostream &os, const sequence &s)
	{
		for(unsigned i=0;i<s.len_;++i)
			os << AlphabetAttributes<_val>::ALPHABET[s.data_[i]];
		return os;
	}
	friend Text_buffer& operator<<(Text_buffer &buf, const sequence &s)
	{
		for(unsigned i=0;i<s.len_;++i)
			buf << AlphabetAttributes<_val>::ALPHABET[s.data_[i]];
		return buf;
	}
	/*friend std::ostream& operator<<(std::ostream &os, const sequence &s)
	{
		std::cout << "co = " << s.clipping_offset_ << std::endl;
		for(unsigned i=s.clipping_offset_;i<s.len_;++i) {
			if(s.data_[i] == 24)
				break;
			os << mask_critical(s.data_[i]);
		}
		return os;
	}*/
	size_t	 len_;
	int		 clipping_offset_;
	_val	*data_;
};


#endif /* SEQUENCE_H_ */
