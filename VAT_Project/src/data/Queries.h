
#ifndef QUERIES_H_
#define QUERIES_H_

#include "../basic/Translator.h"
#include "../tools/complexity_filter.h"
#include "SortedList.h"
#include "../basic/Statistics.h"
#include "SequenceSet.h"
#include <memory>

using std::auto_ptr;

auto_ptr<SeedHistogram> query_hst;
unsigned current_query_chunk;

struct query_source_seqs
{
	static const SequenceSet<DNA>& get()
	{ return *data_; }
	static SequenceSet<DNA> *data_;
};

SequenceSet<DNA>* query_source_seqs::data_ = 0;

template<typename _val>
struct QuerySeqs
{
	static const SequenceSet<_val>& get()
	{ return *data_; }
	static SequenceSet<_val> *data_;
};

template<typename _val> SequenceSet<_val>* QuerySeqs<_val>::data_ = 0;

struct query_ids
{
	static const AlphabetSet<char,0>& get()
	{ return *data_; }
	static AlphabetSet<char,0> *data_;
};

AlphabetSet<char,0>* query_ids::data_ = 0;

#endif /* QUERIES_H_ */
