
#ifndef QUERIES_H_
#define QUERIES_H_

#include "../basic/Translator.h"
#include "../util/complexity_filter.h"
#include "RadixCluster.h"
#include "../basic/statistics.h"
#include "sequence_set.h"

auto_ptr<seed_histogram> query_hst;
unsigned current_query_chunk;

struct query_source_seqs
{
	static const Sequence_set<DNA>& get()
	{ return *data_; }
	static Sequence_set<DNA> *data_;
};

Sequence_set<DNA>* query_source_seqs::data_ = 0;

template<typename _val>
struct query_seqs
{
	static const Sequence_set<_val>& get()
	{ return *data_; }
	static Sequence_set<_val> *data_;
};

template<typename _val> Sequence_set<_val>* query_seqs<_val>::data_ = 0;

struct query_ids
{
	static const String_set<char,0>& get()
	{ return *data_; }
	static String_set<char,0> *data_;
};

String_set<char,0>* query_ids::data_ = 0;

#endif /* QUERIES_H_ */
