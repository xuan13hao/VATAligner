

#ifndef TEXT_BUFFER_H_
#define TEXT_BUFFER_H_

#include <stdlib.h>

struct Text_buffer
{

	Text_buffer():
		data_ (0),
		ptr_ (data_)
	{ }

	void reserve(size_t n)
	{
		const size_t s = ptr_ - data_, new_size = s + n + block_size - ((s+n) & (block_size-1));
		data_ = (char*)realloc(data_, new_size);
		ptr_ = data_ + s;
		if(data_ == 0) throw memory_alloc_exception();
	}

	void operator+=(size_t n)
	{ ptr_ += n; }

	operator char*() const
	{ return ptr_; }

	char* get_begin() const
	{ return data_; }

	void clear()
	{ ptr_ = data_; }

	template<typename _t>
	Text_buffer& write(const _t& data)
	{
		reserve(sizeof(_t));
		*reinterpret_cast<_t*>(ptr_) = data;
		ptr_ += sizeof(_t);
		return *this;
	}

	void write_c_str(const char* s)
	{
		const size_t l = strlen(s)+1;
		reserve(l);
		memcpy(ptr_, s, l);
		ptr_ += l;
	}

	void write_c_str(const char*s, size_t len)
	{
		reserve(len+1);
		memcpy(ptr_, s, len);
		ptr_ += len;
		*(ptr_++) = '\0';
	}

	Text_buffer& write_packed(unsigned x)
	{
		if(x <= (unsigned)std::numeric_limits<uint8_t>::max())
			write((uint8_t)x);
		else if(x <= (unsigned)std::numeric_limits<uint16_t>::max())
			write((uint16_t)x);
		else
			write(x);
		return *this;
	}

	size_t size() const
	{ return ptr_ - data_; }

	Text_buffer& operator<<(const string &s)
	{
		const size_t l = s.length();
		reserve(l);
		memcpy(ptr_, s.c_str(), l);
		ptr_ += l;
		return *this;
	}

	Text_buffer& operator<<(const char* s)
	{
		const size_t l = strlen(s);
		reserve(l);
		memcpy(ptr_, s, l);
		ptr_ += l;
		return *this;
	}

	Text_buffer& operator<<(char c)
	{
		reserve(1);
		*(ptr_++) = c;
		return *this;
	}

	Text_buffer& operator<<(uint32_t x)
	{
		//write(x);
		reserve(16);
		ptr_ += sprintf(ptr_, "%u", x);
		return *this;
	}

	Text_buffer& operator<<(int x)
	{
		//write(x);
		reserve(16);
		ptr_ += sprintf(ptr_, "%i", x);
		return *this;
	}

	Text_buffer& operator<<(double x)
	{
		reserve(32);
		ptr_ += sprintf(ptr_, "%.1lf", x);
		return *this;
	}

	Text_buffer& print_e(double x)
	{
		reserve(32);
		ptr_ += sprintf(ptr_, "%.1le", x);
		return *this;
	}

	/*Text_buffer& operator<<(uint8_t x)
	{
		write(x);
		return *this;
	}*/

	template<typename _t>
	Text_buffer& operator<<(const vector<_t> &v)
	{
		const size_t l = v.size() * sizeof(_t);
		reserve(l);
		memcpy(ptr_, v.data(), l);
		ptr_ += l;
		return *this;
	}

protected:
	enum { block_size = 65536 };
	char *data_, *ptr_;

};

#endif /* TEXT_BUFFER_H_ */
