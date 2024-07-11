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

#ifndef __UTILSTOOLS_H__
#define __UTILSTOOLS_H__

#include <vector>
using std::vector;
using std::string;
using std::pair;

template<typename _it, typename _key>
class Map
{
	public:
	struct Iterator
	{
		Iterator(const _it& begin, const _it& parent_end):
			begin_ (begin),
			parent_end_ (parent_end),
			end_ (get_end())
		{ }
		void operator++()
		{ begin_ = end_; end_ = get_end(); }
		bool valid() const
		{ return begin_ < parent_end_; }
		_it& begin()
		{ return begin_; }
		_it& end()
		{ return end_; }
	private:
		_it get_end() const
		{
			_it i = begin_;
			while(i < parent_end_ && _key()(*i) == _key()(*begin_)) ++i;
			return i;
		}
		_it begin_, parent_end_, end_;
	};

	Map(const _it &begin, const _it &end):
		begin_ (begin),
		end_ (end)
	{ }

	Iterator begin()
	{ return Iterator(begin_, end_); }

private:
	_it begin_, end_;

};

template<typename _t>
struct Double_buffer
{

	inline void init(size_t size, size_t padding, size_t padding_front, _t init)
	{
		const size_t total = size + padding + padding_front;
		data_.resize(total*2);
		ptr1 = &data_[padding_front];
		ptr2 = &data_[total+padding_front];
		for(size_t i=0;i<total*2;++i)
			data_[i] = init;
	}

	inline pair<_t*,_t*> get(int)
	{ std::swap(ptr1, ptr2); return pair<_t*,_t*> (ptr2, ptr1); }

	inline _t* last()
	{ return ptr1; }

private:
	_t *ptr1, *ptr2;
	vector<_t> data_;

};


struct Right { };
struct Left { };

template<typename _val>
_val get_dir(const _val* x, int i, const Right&)
{ return *(x+i); }

template<typename _val>
_val get_dir(const _val* x, int i, const Left&)
{ return *(x-i); }

template<typename _val>
const _val* get_dir_ptr(const _val* x, int i, const Right&)
{ return x+i; }

template<typename _val>
const _val* get_dir_ptr(const _val* x, int i, const Left&)
{ return x-i; }

template<typename _val>
const _val* inc_dir(const _val* x, const Right&)
{ return x+1; }

template<typename _val>
const _val* inc_dir(const _val* x, const Left&)
{ return x-1; }



struct Binary_buffer : public vector<char>
{

	struct Iterator
	{
		Iterator(vector<char>::const_iterator begin, vector<char>::const_iterator end):
			ptr_ (begin),
			end_ (end)
		{ }
		Iterator& operator>>(uint32_t &x)
		{ read(x); return *this; }
		Iterator& operator>>(uint8_t &x)
		{ read(x); return *this; }
		template<typename _t>
		void read(_t &x)
		{
			check(sizeof(_t));
			x = *(_t*)(&*ptr_);
			ptr_ += sizeof(_t);
		}
		template<typename _t>
		void read(vector<_t> &v, size_t count)
		{
			const size_t l = sizeof(_t) * count;
			check(l);
			v.resize(count);
			memcpy(v.data(), &*ptr_, l);
			ptr_ += l;
		}
		void read_packed(uint8_t length, uint32_t &dst)
		{
			switch(length) {
			case 0: uint8_t x; read(x); dst = x; break;
			case 1: uint16_t y; read(y); dst = y; break;
			case 2: read(dst);
			}
		}
		Iterator& operator>>(string &dst)
		{
			dst.clear();
			char c;
			while(read(c), c != '\0')
				dst.push_back(c);
			return *this;
		}
		bool good() const
		{ return ptr_ < end_; }
	private:
		void check(size_t size) const
		{ if(ptr_+size > end_) throw std::runtime_error("Unexpected end of file."); }
		vector<char>::const_iterator ptr_, end_;
	};

	Iterator begin() const
	{ return Iterator (vector<char>::begin(), vector<char>::end()); }

};

class TempFile
{
	public:
	TempFile():
		file_des_ (-1)
	{
		char *name = new char[VATParameters::tmpdir.length()+20];
		sprintf(name, "%s/VAT-tmp-XXXXXX", VATParameters::tmpdir.c_str());
		file_des_ = mkstemp(name);
		if(file_des_ == -1)
			throw std::runtime_error("Error opening temporary file.");
		unlink(name);
		delete[] name;
	}

	int file_descriptor() const
	{ return file_des_; }

private:

	int file_des_;

};



#endif // __UTILSTOOLS_H__