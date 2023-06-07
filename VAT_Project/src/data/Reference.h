

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include <memory>
#include <string>
#include <numeric>
#include "../tools/binary_file.h"
#include "SortedList.h"
#include "../basic/Statistics.h"
#include "../data/SeedHistogram.h"
#include "../tools/hash_table.h"
#include "../tools/hash_function.h"
#include "../basic/PackedLocations.h"
#include "SequenceSet.h"
#include "boost/ptr_container/ptr_vector.hpp"
#include "MaskSequence.h"

using std::auto_ptr;
using boost::ptr_vector;

struct Reference_header
{
	Reference_header():
		unique_id (0x24af8a415ee186dllu),
		build (VATConsts::build_version),
		long_addressing (false),
		sequences (0),
		letters (0)
	{ }
	uint64_t unique_id;
	uint32_t build;
	bool long_addressing;
	unsigned n_blocks;
	size_t sequences, letters;
	double block_size;
#ifdef EXTRA
	Sequence_type sequence_type;
#endif
} ref_header;

struct Database_format_exception : public exception
{
	virtual const char* what() const throw()
	{ return "Database file is not a VAT database."; }
};

template<typename _val>
struct Database_file : public Input_stream
{
	Database_file():
		Input_stream (VATParameters::database)
	{
		if(this->read(&ref_header, 1) != 1)
			throw Database_format_exception ();
		if(ref_header.unique_id != Reference_header ().unique_id)
			throw Database_format_exception ();
		if(ref_header.build > VATConsts::build_version || ref_header.build < VATConsts::build_compatibility)
			throw invalid_database_version_exception();
#ifdef EXTRA
		if(sequence_type(_val()) != ref_header.sequence_type)
			throw std::runtime_error("Database has incorrect sequence type for current alignment mode.");
#endif
	}
	void rewind()
	{ this->seek(sizeof(Reference_header)); }
};

template<typename _val>
struct ReferenceSeqs
{
	static const Masked_sequence_set<_val>& get()
	{ return *(const Masked_sequence_set<_val>*)data_; }
	static Masked_sequence_set<_val>& get_nc()
	{ return *(Masked_sequence_set<_val>*)data_; }
	static SequenceSet<_val> *data_;
};

template<typename _val> SequenceSet<_val>* ReferenceSeqs<_val>::data_ = 0;

struct ReferenceIds
{
	static const AlphabetSet<char,0>& get()
	{ return *data_; }
	static AlphabetSet<char,0> *data_;
};

AlphabetSet<char,0>* ReferenceIds::data_ = 0;

SeedHistogram ref_hst;
unsigned current_ref_block;

size_t max_id_len(const AlphabetSet<char,0> &ids)
{
	size_t max (0);
	for(size_t i=0;i<ids.get_length(); ++i)
		max = std::max(max, find_first_of(ids[i].c_str(), VATConsts::id_delimiters));
	return max;
}

struct Ref_map
{
	Ref_map():
		next_ (0)
	{ }
	void init(unsigned ref_count)
	{
		const unsigned block = current_ref_block;
		if(data_.size() < block+1) {
			data_.resize(block+1);
			data_[block].insert(data_[block].end(), ref_count, std::numeric_limits<unsigned>::max());
		}
	}
	template<typename _val>
	uint32_t get(unsigned block, unsigned i)
	{
		uint32_t n = data_[block][i];
		if(n != std::numeric_limits<unsigned>::max())
			return n;
		else {
			tthread::lock_guard<tthread::mutex> lock (mtx_);
			n = data_[block][i];
			if(n != std::numeric_limits<uint32_t>::max())
				return n;
			n = next_++;
			data_[block][i] = n;
			len_.push_back(ReferenceSeqs<_val>::get().length(i));
			name_.push_back(get_str(ReferenceIds::get()[i].c_str(), VATConsts::id_delimiters));
			return n;
		}
	}
	/*template<typename _val>
	void finish()
	{
		vector<pair<unsigned,unsigned> > v;
		for(unsigned i=0;i<data_.size();++i)
			if(data_[i] != std::numeric_limits<unsigned>::max())
				v.push_back(pair<unsigned,unsigned> (data_[i], i));
		std::sort(v.begin(), v.end());
		for(vector<pair<unsigned,unsigned> >::const_iterator i = v.begin(); i!=v.end(); ++i) {
			const char* s = ref_ids::get()[i->second].c_str();
			buf_ << (uint32_t)(strlen(s)+1);
			buf_.write_c_str(s);
			buf_ << (uint32_t)ref_seqs<_val>::get()[i->second].length();
		}
	}*/
private:
	tthread::mutex mtx_;
	vector<vector<uint32_t> > data_;
	vector<uint32_t> len_;
	ptr_vector<string> name_;
	uint32_t next_;
	friend struct VATOutput;
} ref_map;

#endif /* REFERENCE_H_ */
