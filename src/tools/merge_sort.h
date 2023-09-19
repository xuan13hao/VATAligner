
#ifndef MERGE_SORT_H_
#define MERGE_SORT_H_

#include <algorithm>
#include "thread.h"

template<typename _it>
void merge_sort(_it begin, _it end, unsigned n_threads, unsigned level = 0)
{
	ptrdiff_t diff = end - begin;
	if(diff <= 1)
		return;

	if(1u << level >= n_threads) {
		std::sort(begin, end);
		return;
	}

	_it mid = begin + diff/2;
	thread *left = launch_thread(merge_sort<_it>, begin, mid, n_threads, level+1);
	thread *right = launch_thread(merge_sort<_it>, mid, end, n_threads, level+1);
	left->join();
	right->join();
	delete left;
	delete right;
	std::inplace_merge(begin, mid, end);
}

#endif /* MERGE_SORT_H_ */
