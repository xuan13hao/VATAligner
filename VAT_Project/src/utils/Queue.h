#ifndef __QUEUE_H__
#define __QUEUE_H__



#include <limits>
#include <mutex>
#include <condition_variable>

class Queue
{
	public:
	enum { end = size_t(-1) };
	Queue(size_t begin, size_t end) :
		next_(begin),
		block_(false),
		end_(end)
	{}
	template<typename _f>
	size_t get(_f &f)
	{
		std::unique_lock<std::mutex> lock(mtx_);
		while (block_)
			cond_.wait(lock);
		const size_t q = next_++;
		if (q >= end_) {
			return Queue::end;
		}
		block_ = f(q);
		return q;
	}
	size_t next() const
	{
		return next_;
	}
	size_t get_end() const
	{
		return end_;
	}
	void release() {
		block_ = false;
		cond_.notify_all();
	}
private:
	std::mutex mtx_;
	std::condition_variable cond_;
	volatile size_t next_;
	volatile bool block_;
	const size_t end_;
};
#endif // __QUEUE_H__