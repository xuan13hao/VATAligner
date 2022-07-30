
#ifndef INTERVAL_H_
#define INTERVAL_H_

#include <iostream>
#include <algorithm>

struct interval
{
	interval() :
		begin_(0),
		end_(0)
	{ }
	interval(int begin, int end) :
		begin_(begin),
		end_(end)
	{ }
	int length() const
	{
		return end_ > begin_ ? end_ - begin_ : 0;
	}
	unsigned overlap(const interval &rhs) const
	{
		return intersect(*this, rhs).length();
	}
	double overlap_factor(const interval &rhs) const
	{
		return (double)overlap(rhs) / (double)length();
	}
	bool includes(int p) const
	{
		return p >= begin_ && p < end_;
	}
	bool contains(const interval &i) const
	{
		return begin_ <= i.begin_ && end_ >= i.end_;
	}
	friend std::ostream& operator<<(std::ostream &os, const interval &x)
	{
		os << "[" << x.begin_ << ";" << x.end_ << "]"; return os;
	}
	bool operator<(const interval &rhs) const
	{
		return begin_ < rhs.begin_;
	}
	void merge(const interval &k)
	{
		begin_ = std::min(begin_, k.begin_);
		end_ = std::max(end_, k.end_);
	}
	friend interval intersect(const interval &lhs, const interval &rhs);
	int begin_, end_;
};

inline interval intersect(const interval &lhs, const interval &rhs)
{
	return interval(std::max(lhs.begin_, rhs.begin_), std::min(lhs.end_, rhs.end_));
}

#endif