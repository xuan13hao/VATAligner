

#ifndef STRING_SET_H_
#define STRING_SET_H_

#include <vector>

using std::vector;

template<typename _t, char _pchar = 0xff, size_t _padding = 1>
class AlphabetSet
{

	public:
	static const unsigned PERIMETER_PADDING = 1;
	static const _t PADDING_CHAR;

	AlphabetSet():
		data_ (PERIMETER_PADDING)
	{ limits_.push_back(PERIMETER_PADDING); }

	void finish_reserve()
	{
		data_.resize(raw_len() + PERIMETER_PADDING);
		for(unsigned i=0;i<PERIMETER_PADDING;++i) {
			data_[i] = PADDING_CHAR;
			data_[raw_len()+i] = PADDING_CHAR;
		}
	}

	void push_back(const vector<_t> &v)
	{
		limits_.push_back(raw_len() + v.size() + _padding);
		data_.insert(data_.end(), v.begin(), v.end());
		data_.insert(data_.end(), _padding, PADDING_CHAR);
	}

	void fill(size_t n, _t v)
	{
		limits_.push_back(raw_len() + n + _padding);
		data_.insert(data_.end(), n, v);
		data_.insert(data_.end(), _padding, PADDING_CHAR);
	}

	_t* ptr(size_t i)
	{ return &data_[limits_[i]]; }

	const _t* ptr(size_t i) const
	{ return &data_[limits_[i]]; }

	size_t length(size_t i) const
	{ return limits_[i+1] - limits_[i] - _padding; }

	size_t get_length() const
	{ return limits_.size() - 1; }

	void save(OutputStreamer &file) const
	{
		file.write(limits_);
		file.write(data_);
	}

	AlphabetSet(Input_stream &file)
	{
		file.read(limits_);
		file.read(data_);
	}

	static void skip(Input_stream &file)
	{
		file.skip_vector<size_t>();
		file.skip_vector<_t>();
	}
	//raw data size
	size_t raw_len() const
	{ return limits_.back(); }

	size_t letters() const
	{ return raw_len() - get_length() - PERIMETER_PADDING; }

	_t* data(ptrdiff_t p = 0)
	{ return &data_[p]; }

	const _t* data(ptrdiff_t p = 0) const
	{ return &data_[p]; }

	size_t position(const _t* p) const
	{ return p - data(); }

	size_t position(size_t i, size_t j) const
	{ return limits_[i] + j; }

	std::pair<size_t,size_t> local_position(size_t p) const
	{
		//i is the query number 
		size_t i = std::upper_bound(limits_.begin(), limits_.end(), p) - limits_.begin() - 1;
		return std::pair<size_t,size_t> (i, p - limits_[i]);
	}

	sequence<const _t> operator[](size_t i) const
	{ return sequence<const _t> (ptr(i), length(i)); }

	sequence<_t> operator[](size_t i)
	{ return sequence<_t> (ptr(i), length(i)); }

private:

	vector<_t> data_;
	vector<size_t> limits_;

};

template<typename _t, char _pchar, size_t _padding> const _t AlphabetSet<_t,_pchar,_padding>::PADDING_CHAR = _pchar;

#endif /* STRING_SET_H_ */
