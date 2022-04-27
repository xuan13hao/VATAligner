
#ifndef QUERIES_H_
#define QUERIES_H_

#include "../basic/Translator.h"
#include "../util/complexity_filter.h"
#include "RadixCluster.h"
#include "../basic/statistics.h"
#include "SequenceSet.h"

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
struct query_seqs
{
	static const SequenceSet<_val>& get()
	{ return *data_; }
	static SequenceSet<_val> *data_;
};

template<typename _val> SequenceSet<_val>* query_seqs<_val>::data_ = 0;

struct query_ids
{
	static const AlphabetSet<char,0>& get()
	{ return *data_; }
	static AlphabetSet<char,0> *data_;
};

AlphabetSet<char,0>* query_ids::data_ = 0;

#endif /* QUERIES_H_ */
