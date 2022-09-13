#ifndef __FLAT_ARRAY_H__
#define __FLAT_ARRAY_H__

#include <vector>

template<typename _t>
class FlatArray {
	public:
	FlatArray() {
		limits_.push_back(0);
	}

	void push_back(const _t &x) {
		data_.push_back(x);
		++limits_.back();
	}

	void next() {
		limits_.push_back(limits_.back());
	}

	void clear() {
		data_.clear();
		limits_.clear();
		limits_.push_back(0);
	}

	size_t size() const {
		return limits_.size() - 1;
	}

	const _t* begin(size_t i) const {
		return &data_[limits_[i]];
	}

	const _t* end(size_t i) const {
		return &data_[limits_[i + 1]];
	}

private:

	std::vector<_t> data_;
	std::vector<size_t> limits_;

};

#endif // __FLAT_ARRAY_H__