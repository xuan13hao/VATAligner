#ifndef SORTED_LIST_H_
#define SORTED_LIST_H_

#include <memory>
#include <map>
#include "../tools/sort_simd/include/avx256/simd_sort.h"
#include "../tools/util.h"
#include "SeedHistogram.h"
#include "../Commons/PackedLocations.h"
#define MAX_uint64 9999999

using std::auto_ptr;
using namespace std;

std::pair<int, std::string> compute_minimizer(const std::string& kmer, int w) {
    int k = kmer.length();
    
    if (k <= w) {
        // If k is less than or equal to w, the entire k-mer is the minimizer.
        return {0, kmer};
    }

    std::string minimizer = kmer.substr(0, w);  // Initialize with the first w characters
    int minPosition = 0;

    for (int i = w; i < k; ++i) {
        std::string currentKmer = kmer.substr(i - w + 1, w);  // Current w-mer within the window

        // Update minimizer if the current w-mer is lexicographically smaller
        if (currentKmer < minimizer) {
            minimizer = currentKmer;
            minPosition = i - w + 1;  // Update the position of the minimizer
        }
    }

    return {minPosition, minimizer};
}
 	// Function to compute the hash value for a k-mer
	int computeKmerHash(const std::string& kmer) {
		int hash = 0;

		for (char nucleotide : kmer) {
			if (nucleotide == 'A') {
				hash = hash * 10 + 1;
			} else if (nucleotide == 'C') {
				hash = hash * 10 + 2;
			} else if (nucleotide == 'G') {
				hash = hash * 10 + 3;
			} else if (nucleotide == 'T') {
				hash = hash * 10 + 4;
			}
		}

		return hash;
	}
std::map<int, std::string> generateKmerMap(const std::string& sequence, int k) {
    std::map<int, std::string> kmerMap; // Key is the sum, and value is the k-mer
    int sequenceLength = sequence.length();

    // The total number of possible k-mers is 4^k
    int totalKmers = 1;
    for (int i = 0; i < k; ++i) {
        totalKmers *= 4;
    }

    for (int i = 0; i < totalKmers; ++i) {
        std::string kmer;

        // Generate a k-mer based on the current index 'i'
        int index = i;
        int sum = 0; // Initialize the sum for this k-mer

        for (int j = 0; j < k; ++j) {
            int nucleotideIndex = index % 4; // There are 4 nucleotides (A, C, G, T)
            char nucleotide;

            // Map nucleotides to numbers
            if (nucleotideIndex == 0) {
                nucleotide = 'A';
                sum = sum * 10 + 1;
            } else if (nucleotideIndex == 1) {
                nucleotide = 'C';
                sum = sum * 10 + 2;
            } else if (nucleotideIndex == 2) {
                nucleotide = 'G';
                sum = sum * 10 + 3;
            } else {
                nucleotide = 'T';
                sum = sum * 10 + 4;
            }

            kmer += nucleotide; // Build the k-mer
            index /= 4;
        }

        // Store the k-mer as a value, and the sum as the key
        kmerMap[sum] = kmer;
    }

    return kmerMap;
}

template<typename _val>    
std::string to_string(const sequence<const _val> &s) 
{
	std::stringstream ss;
	for(unsigned i=0;i<s.len_;++i)
	{
		ss << AlphabetAttributes<_val>::ALPHABET[s.data_[i]];
	}
	return ss.str();
}

template<typename _pos>
class SortedList
{
	public:
	typedef SortedList<typename packed_sequence_location<_pos>::type> Type;

	class Tuple
	{
		public:
		Tuple():
			key (),
			value ()
		{ 

		}
		Tuple(unsigned key, _pos value):
			key (key),
			value (value)
		{ 

		}
		Tuple(unsigned key, _pos value,unsigned suffix_bit):
			key (key),
			value (value),
			suffix_bit (suffix_bit)
		{ 

		}
		void set_key(unsigned i)
		{
			this->key = i;
		}
		void set_suffix(unsigned i)
		{
			this->suffix_bit = i;
		}
		void set_value(_pos i)
		{
			this->value = i;
		}
		unsigned get_key()
		{
			unsigned k = this->key;
			return k;
		}
		unsigned get_suffix()
		{
			unsigned k = this->suffix_bit;
			return k;
		}
		bool operator<(const Tuple &rhs) const
		{ 
			return key < rhs.key; 
		}
		// entry operator++() const
		// { return key < rhs.key; }
		unsigned	key;
		unsigned	suffix_bit;
		_pos		value;
	} __attribute__((packed));////Tells the compiler to allocate the variable x at a 16-byte aligned memory address instead of the default 4-byte alignment.

	static char* alloc_buffer(const SeedHistogram &hst)
	{ 
		return new char[sizeof(Tuple) * hst.max_chunk_size()]; 
	}



	template<typename _val>
	SortedList(char *buffer, const SequenceSet<_val> &seqs, const Shape &sh, const ShapeHistogram &hst, const seedp_range &range):
		limits_ (hst, range),
		data_ (reinterpret_cast<Tuple*>(buffer))
	{
		TimerTools timer ("Building seed list", false);
		Build_context<_val> build_context (seqs, sh, range, build_iterators(hst));
		launch_scheduled_thread_pool(build_context, VATConsts::seqp, 2*VATParameters::threads());
		timer.go("Sorting seed list");
		Sort_context sort_context (*this);
		launch_scheduled_thread_pool(sort_context, VATConsts::seedp, 2*VATParameters::threads());
	}

	template<typename _t>
	struct Iterator_base
	{
		Iterator_base(_t *i, _t *end):
			i (i),
			end (end),
			n (count())
		{ }
		size_t count() const
		{
			_t *k (i);
			size_t n (0);
			while(k < end && (k++)->key == i->key)
				++n;
			return n;
		}
		void operator++()
		{ i += n; n = count(); }
		_pos operator[](unsigned k) const
		{ return (i+k)->value; }
		_t* get(unsigned k)
		{ return i+k; }
		bool at_end() const
		{ return i >= end; }
		unsigned key() const
		{ return i->key; }
		_t *i, *end;
		size_t n;
	};

	typedef Iterator_base<Tuple> iterator;
	typedef Iterator_base<const Tuple> const_iterator;

	const_iterator get_partition_cbegin(unsigned p) const
	{ return const_iterator (cptr_begin(p), cptr_end(p)); }

	iterator get_partition_begin(unsigned p) const
	{ return iterator (ptr_begin(p), ptr_end(p)); }

private:

	typedef Static_matrix<Tuple*,VATConsts::seqp,VATConsts::seedp> Ptr_set;

	struct buffered_iterator
	{
		static const unsigned BUFFER_SIZE = 16;
		buffered_iterator(Tuple **ptr)
		{
			memset(n, 0, sizeof(n));
			memcpy(this->ptr, ptr, sizeof(this->ptr));
		}
		//insert seed with postion into seedp_range
		void push(seed key, _pos value, const seedp_range &range)
		{
			const unsigned p (seed_partition(key)); //63505
			if(range.contains(p)) {
				assert(n[p] < BUFFER_SIZE);
				// buf[p][n[p]++] = Tuple (key, value);
				buf[p][n[p]++] = Tuple (seed_partition_offset(key), value);
				if(n[p] == BUFFER_SIZE)
					flush(p);
			}
			// const unsigned p (seed_partition(key));
			// if (range.contains(p)) {
			// 	assert(n[p] < BUFFER_SIZE);
				
			// 	// Check the distance between the keys of the current Tuple and the new Tuple 63501
			// 	if (n[p] > 0 && key - buf[p][n[p] - 1].key > 5) {
			// 		// If the distance is smaller than 5, update the value of the existing Tuple
			// 		buf[p][n[p] - 1].value = value;
			// 	} else {
			// 		// Otherwise, add the new Tuple to the buffer
			// 		buf[p][n[p]++] = Tuple(seed_partition_offset(key), value);
			// 	}
				
			// 	if (n[p] == BUFFER_SIZE)
			// 		flush(p);
			// }

		}
		void flush(unsigned p)
		{
			//when buf size > 16, copy buf to ptr
			memcpy(ptr[p], buf[p], n[p] * sizeof(Tuple));
			ptr[p] += n[p];//mv point with n[p] bits
			n[p] = 0;
		}
		void flush()
		{
			for(unsigned p=0;p<VATConsts::seedp;++p)
				if(n[p] > 0)
					flush(p);
		}
		Tuple* ptr[VATConsts::seedp];
		Tuple  	 buf[VATConsts::seedp][BUFFER_SIZE];
		uint8_t  n[VATConsts::seedp];
	};

	Tuple* ptr_begin(unsigned i) const
	{ 
		// cout<<"data_[limits_[i]] = "<<&data_[limits_[i]]<<endl;
		// for (size_t j = 0; j < limits_[i].size(); j++)
		// {
		// 	cout<<"limit = "<<limits_[i][j]<<endl;
		// }
		
		return &data_[limits_[i]]; 
	}

	Tuple* ptr_end(unsigned i) const
	{ return &data_[limits_[i+1]]; }

	const Tuple* cptr_begin(unsigned i) const
	{ return &data_[limits_[i]]; }

	const Tuple* cptr_end(unsigned i) const
	{ return &data_[limits_[i+1]]; }

	template<typename _val>
	struct Build_context
	{
		Build_context(const SequenceSet<_val> &seqs, const Shape &sh, const seedp_range &range, Ptr_set *iterators):
			seqs (seqs),
			sh (sh),
			range (range),
			iterators (iterators),
			seq_partition (seqs.partition())
		{ }
		void operator()(unsigned thread_id, unsigned seqp) const
		{
			build_seqp<_val>(seqs,
					seq_partition[seqp],
					seq_partition[seqp+1],
					(*iterators)[seqp],
					sh,
					range);
		}
		const SequenceSet<_val> &seqs;
		const Shape &sh;
		const seedp_range &range;
		const auto_ptr<Ptr_set> iterators;
		const vector<size_t> seq_partition;
	};

	template<typename _val>
	static void build_seqp(const SequenceSet<_val> &seqs, unsigned begin, unsigned end, Tuple **ptr, const Shape &sh, const seedp_range &range)
	{
		uint64_t key;
		//init buffered iterator via entry size
		// cout<<"sh.length_+1 = "<<sh.length_+1<<endl;
		auto_ptr<buffered_iterator> it (new buffered_iterator(ptr));
		for(size_t i=begin;i<end;++i) 
		{
			const sequence<const _val> seq = seqs[i];
			std::string str_seq = to_string(seq);
			int w = VATParameters::window;
			int k = sh.length_;
			// for(unsigned j=0;j<seq.length()-sh.length_+1; ++j) 
			// {
			// 	std::string k_mer = str_seq.substr(j, k);
			// 	auto result = compute_minimizer(k_mer,w);
			// 	key = computeKmerHash(result.second);
			// 	unsigned pos = j + result.first;
			// 	it->push(key, seqs.position(i, pos), range);
			// }
			if(seq.length()<sh.length_) continue;
	
			for(unsigned j=0;j<seq.length()-sh.length_+1; ++j) 
			{
				// std::string window = str_seq.substr(j, w);
				// pre_kmer = computeKmerHash(window);
				if(sh.set_seed(key, &seq[j]))//get key via seq
					it->push(key, seqs.position(i, j), range);
					
			}
		}
		//copy all buffer data into ptr set
		it->flush();
	}

	Ptr_set* build_iterators(const ShapeHistogram &hst) const
	{
		Ptr_set *iterators = new Ptr_set;
		for(unsigned i=0;i<VATConsts::seedp;++i)
			(*iterators)[0][i] = ptr_begin(i);

		for(unsigned i=1;i<VATConsts::seqp;++i) {
			for(unsigned j=0;j<VATConsts::seedp;++j)
			{
				(*iterators)[i][j] = (*iterators)[i-1][j] + hst[i-1][j];
				//cout<<"iterator = "<<(*iterators)[i][j]<<", hst = "<<(hst[i-1][j])<<endl;
			}
				

		}
		return iterators;
	}
	void sortSIMD(Tuple* s ,Tuple* e ) 
	// void sortSIMD() 
	{
		uint64_t N = e - s;
		uint64_t p = N % 8;
		uint64_t n = N;
		std::pair<uint64_t, uint64_t> *rand_arr;
		std::pair<uint64_t, uint64_t> *soln_arr;
		if(p != 0)
		{   
			n = N + 8 - p  ;
		}
		// allocate n memory
		aligned_init<std::pair<uint64_t, uint64_t> >(rand_arr, n);
		for (size_t i = 0; i < n; i++)
		{
			if (i < N )
			{
				uint64_t k = (uint64_t)s[i].key;
				uint64_t v = (uint64_t)s[i].value;
				// std::pair<uint64_t ,uint64_t> pair_(k,v);
				rand_arr[i].first = k;
				rand_arr[i].second = v;
			}else
			{
				// std::pair<uint64_t ,uint64_t> Padding(MAX_uint64,MAX_uint64);
				rand_arr[i].first = MAX_uint64;
				rand_arr[i].second = MAX_uint64;
			}	
		}
   		aligned_init<std::pair<uint64_t ,uint64_t> >(soln_arr, n);
		std::copy(rand_arr, rand_arr + n, soln_arr);
		avx2::SIMDSort(n, soln_arr);
		for (size_t i = 0; i < n; i++)
    	{
			if ((soln_arr[i].first != MAX_uint64) && (soln_arr[i].second != MAX_uint64))
			{
				unsigned k = (unsigned)soln_arr[i].first;
				_pos v = (_pos)soln_arr[i].second;
				s[i].key = k;
				s[i].value = v;
			}
    	}
		delete rand_arr;
		delete soln_arr;
	}

	struct Sort_context
	{
		Sort_context(SortedList &sl):
			sl (sl)
		{ }

		bool is_power_of_two(int n) {
			return (n != 0) && ((n & (n - 1)) == 0);
		}
		void operator()(unsigned thread_id ,unsigned seedp) const
		{
			int n = sl.ptr_end(seedp) - sl.ptr_begin(seedp) ;
			if((n >= 8 && ((n != 0) && ((n & (n - 1)) == 0))) && VATParameters::simd_sort)
			{
				// cout<<"Using SIMD Sort "<<n<<endl;
				sl.sortSIMD(sl.ptr_begin(seedp), sl.ptr_end(seedp));
			}else
			{
				
				std::sort(sl.ptr_begin(seedp), sl.ptr_end(seedp)); 
				
			}

		}
		SortedList &sl;
	};

	struct Limits : vector<size_t>
	{
		Limits(const ShapeHistogram &hst, const seedp_range &range)
		{
			TimerTools timer ("Computing limits", false);
			this->push_back(0);
			for(unsigned i=0;i<VATConsts::seedp;++i) {
#ifdef EXTRA
				log_stream << i << ' ' << partition_size(hst, i) << endl;
#endif
				this->push_back(this->operator[](i) + (range.contains(i) ? partition_size(hst, i) : 0));
				// cout<<"limit = "<<this->operator[](i)<<endl;
			}
		}
	};

	const Limits limits_;
	Tuple *data_;
};

#endif /* SORTED_LIST_H_ */

// #ifndef SORTED_LIST_H_
// #define SORTED_LIST_H_

// #include <memory>
// #include <map>
// #include "../tools/sort_simd/include/avx256/simd_sort.h"
// #include "../tools/util.h"
// #include "SeedHistogram.h"
// #include "../Commons/PackedLocations.h"
// #define MAX_uint64 9999999

// using std::auto_ptr;
// using namespace std;

// std::pair<int, std::string> compute_minimizer(const std::string& kmer, int w) {
//     int k = kmer.length();
    
//     if (k <= w) {
//         // If k is less than or equal to w, the entire k-mer is the minimizer.
//         return {0, kmer};
//     }

//     std::string minimizer = kmer.substr(0, w);  // Initialize with the first w characters
//     int minPosition = 0;

//     for (int i = w; i < k; ++i) {
//         std::string currentKmer = kmer.substr(i - w + 1, w);  // Current w-mer within the window

//         // Update minimizer if the current w-mer is lexicographically smaller
//         if (currentKmer < minimizer) {
//             minimizer = currentKmer;
//             minPosition = i - w + 1;  // Update the position of the minimizer
//         }
//     }

//     return {minPosition, minimizer};
// }
//  	// Function to compute the hash value for a k-mer
// 	int computeKmerHash(const std::string& kmer) {
// 		int hash = 0;

// 		for (char nucleotide : kmer) {
// 			if (nucleotide == 'A') {
// 				hash = hash * 10 + 1;
// 			} else if (nucleotide == 'C') {
// 				hash = hash * 10 + 2;
// 			} else if (nucleotide == 'G') {
// 				hash = hash * 10 + 3;
// 			} else if (nucleotide == 'T') {
// 				hash = hash * 10 + 4;
// 			}
// 		}

// 		return hash;
// 	}
// std::map<int, std::string> generateKmerMap(const std::string& sequence, int k) {
//     std::map<int, std::string> kmerMap; // Key is the sum, and value is the k-mer
//     int sequenceLength = sequence.length();

//     // The total number of possible k-mers is 4^k
//     int totalKmers = 1;
//     for (int i = 0; i < k; ++i) {
//         totalKmers *= 4;
//     }

//     for (int i = 0; i < totalKmers; ++i) {
//         std::string kmer;

//         // Generate a k-mer based on the current index 'i'
//         int index = i;
//         int sum = 0; // Initialize the sum for this k-mer

//         for (int j = 0; j < k; ++j) {
//             int nucleotideIndex = index % 4; // There are 4 nucleotides (A, C, G, T)
//             char nucleotide;

//             // Map nucleotides to numbers
//             if (nucleotideIndex == 0) {
//                 nucleotide = 'A';
//                 sum = sum * 10 + 1;
//             } else if (nucleotideIndex == 1) {
//                 nucleotide = 'C';
//                 sum = sum * 10 + 2;
//             } else if (nucleotideIndex == 2) {
//                 nucleotide = 'G';
//                 sum = sum * 10 + 3;
//             } else {
//                 nucleotide = 'T';
//                 sum = sum * 10 + 4;
//             }

//             kmer += nucleotide; // Build the k-mer
//             index /= 4;
//         }

//         // Store the k-mer as a value, and the sum as the key
//         kmerMap[sum] = kmer;
//     }

//     return kmerMap;
// }

// 	template<typename _val>    
// 	std::string to_string(const sequence<const _val> &s) 
// 	{
// 		std::stringstream ss;
// 		for(unsigned i=0;i<s.len_;++i)
// 		{
// 			ss << AlphabetAttributes<_val>::ALPHABET[s.data_[i]];
// 		}
// 		return ss.str();
// 	}

// template<typename _pos>
// class SortedTuples
// {
// 	public:
// 	typedef SortedTuples<typename packed_sequence_location<_pos>::type> Type;

// 	class Tuple
// 	{
// 		public:
// 		Tuple():
// 			key (),
// 			value ()
// 		{ 

// 		}
// 		Tuple(unsigned key, _pos value):
// 			key (key),
// 			value (value)
// 		{ 

// 		}
// 		// Tuple(unsigned key, _pos value, unsigned pre_kmer):
// 		// key (key),
// 		// value (value),
// 		// pre_kmer(pre_kmer)
// 		// { 

// 		// }
// 		void set_key(unsigned i)
// 		{
// 			this->key = i;
// 		}
// 		// void set_prekmer(unsigned i)
// 		// {
// 		// 	this->pre_kmer = i;
// 		// }
// 		void set_value(_pos i)
// 		{
// 			this->value = i;
// 		}
// 		unsigned get_key()
// 		{
// 			unsigned k = this->key;
// 			return k;
// 		}
// 		// unsigned get_prekmer()
// 		// {
// 		// 	unsigned k = this->pre_kmer;
// 		// 	return k;
// 		// }
// 		bool operator<(const Tuple &rhs) const
// 		{ 
// 			return key < rhs.key; 
// 		}
// 		// entry operator++() const
// 		// { return key < rhs.key; }
// 		// unsigned	pre_kmer;
// 		unsigned	key;
// 		_pos		value;
// 	} __attribute__((packed));////Tells the compiler to allocate the variable x at a 16-byte aligned memory address instead of the default 4-byte alignment.

// 	static char* alloc_buffer(const SeedHistogram &hst)
// 	{ 
// 		return new char[sizeof(Tuple) * hst.max_chunk_size()]; 
// 	}



// 	template<typename _val>
// 	SortedTuples(char *buffer, const SequenceSet<_val> &seqs, const Shape &sh, const ShapeHistogram &hst, const seedp_range &range):
// 		limits_ (hst, range),
// 		data_ (reinterpret_cast<Tuple*>(buffer))
// 	{
// 		TimerTools timer ("Building seed list", false);
// 		Build_context<_val> build_context (seqs, sh, range, build_iterators(hst));
// 		//2*VATParameters::threads()
// 		launch_scheduled_thread_pool(build_context, VATConsts::seqp, 2*VATParameters::threads());
// 		timer.go("Sorting seed list");
// 		Sort_context sort_context (*this);
// 		launch_scheduled_thread_pool(sort_context, VATConsts::seedp, 2*VATParameters::threads());
// 	}

// 	template<typename _t>
// 	struct Iterator_base
// 	{
// 		Iterator_base(_t *i, _t *end):
// 			i (i),
// 			end (end),
// 			n (count())
// 		{ }
// 		size_t count() const
// 		{
// 			_t *k (i);
// 			size_t n (0);
// 			while(k < end && (k++)->key == i->key)
// 				++n;
// 			return n;
// 		}
// 		void operator++()
// 		{ i += n; n = count(); }
// 		_pos operator[](unsigned k) const
// 		{ return (i+k)->value; }
// 		_t* get(unsigned k)
// 		{ return i+k; }
// 		bool at_end() const
// 		{ return i >= end; }
// 		unsigned key() const
// 		{ return i->key; }
// 		unsigned pre_kmer() const
// 		{ return i->pre_kmer; }
// 		_t *i, *end;
// 		size_t n;
// 	};

// 	typedef Iterator_base<Tuple> iterator;
// 	typedef Iterator_base<const Tuple> const_iterator;

// 	const_iterator get_partition_cbegin(unsigned p) const
// 	{ return const_iterator (cptr_begin(p), cptr_end(p)); }

// 	iterator get_partition_begin(unsigned p) const
// 	{ return iterator (ptr_begin(p), ptr_end(p)); }

// private:

// 	typedef Static_matrix<Tuple*,VATConsts::seqp,VATConsts::seedp> Ptr_set;

// 	struct buffered_iterator
// 	{
// 		static const unsigned BUFFER_SIZE = 16;
// 		buffered_iterator(Tuple **ptr)
// 		{
// 			memset(n, 0, sizeof(n));
// 			memcpy(this->ptr, ptr, sizeof(this->ptr));
// 		}
// 		//insert seed with postion into seedp_range
// 		void push(seed key, _pos value, const seedp_range &range)
// 		{
// 			const unsigned p (seed_partition(key)); //63505
// 			if(range.contains(p)) {
// 				assert(n[p] < BUFFER_SIZE);
// 				buf[p][n[p]++] = Tuple (seed_partition_offset(key), value);
// 				// buf[p][n[p]++] = Tuple (seed_partition_offset(key),value,pre_kmer);
// 				if(n[p] == BUFFER_SIZE)
// 					flush(p);
// 			}

// 		}
// 		void flush(unsigned p)
// 		{
// 			//when buf size > 16, copy buf to ptr
// 			memcpy(ptr[p], buf[p], n[p] * sizeof(Tuple));
// 			ptr[p] += n[p];//mv point with n[p] bits
// 			n[p] = 0;
// 		}
// 		void flush()
// 		{
// 			for(unsigned p=0;p<VATConsts::seedp;++p)
// 				if(n[p] > 0)
// 					flush(p);
// 		}
// 		Tuple* ptr[VATConsts::seedp];
// 		Tuple  	 buf[VATConsts::seedp][BUFFER_SIZE];
// 		uint8_t  n[VATConsts::seedp];
// 	};

// 	Tuple* ptr_begin(unsigned i) const
// 	{ 
// 		return &data_[limits_[i]]; 
// 	}

// 	Tuple* ptr_end(unsigned i) const
// 	{ return &data_[limits_[i+1]]; }

// 	const Tuple* cptr_begin(unsigned i) const
// 	{ return &data_[limits_[i]]; }

// 	const Tuple* cptr_end(unsigned i) const
// 	{ return &data_[limits_[i+1]]; }

// 	template<typename _val>
// 	struct Build_context
// 	{
// 		Build_context(const SequenceSet<_val> &seqs, const Shape &sh, const seedp_range &range, Ptr_set *iterators):
// 			seqs (seqs),
// 			sh (sh),
// 			range (range),
// 			iterators (iterators),
// 			seq_partition (seqs.partition())
// 		{ }
// 		void operator()(unsigned thread_id, unsigned seqp) const
// 		{
// 			build_seqp<_val>(seqs,
// 					seq_partition[seqp],
// 					seq_partition[seqp+1],
// 					(*iterators)[seqp],
// 					sh,
// 					range);
// 		}
// 		const SequenceSet<_val> &seqs;
// 		const Shape &sh;
// 		const seedp_range &range;
// 		const auto_ptr<Ptr_set> iterators;
// 		const vector<size_t> seq_partition;
// 	};

// 	template<typename _val>
// 	static void build_seqp(const SequenceSet<_val> &seqs, unsigned begin, unsigned end, Tuple **ptr, const Shape &sh, const seedp_range &range)
// 	{
// 		uint64_t key;
// 		int w = 3;
// 		//init buffered iterator via entry size
// 		auto_ptr<buffered_iterator> it (new buffered_iterator(ptr));
// 		for(size_t i=begin;i<end;++i) 
// 		{
// 			const sequence<const _val> seq = seqs[i];
// 			std::string str_seq = to_string(seq);
// 			int w = 10;
// 			int k = sh.length_;
// 			// // std::pair<int, std::string> minimizer = findMinimizer()
// 			// for(unsigned j=0;j<seq.length()-sh.length_+1; ++j) 
// 			// {
// 			// 	std::string k_mer = str_seq.substr(j, k);
// 			// 	auto result = compute_minimizer(k_mer,w);
// 			// 	key = computeKmerHash(result.second);
// 			// 	unsigned pos = j + result.first;
// 			// 	it->push(key, seqs.position(i, pos), range);
// 			// }
// 			if(seq.length()<sh.length_) continue;
// 			uint64_t prev_key = 0;
// 			unsigned pre_kmer = 0;
// 			// cout<<"str_seq = "<<str_seq<<endl;
// 			for(unsigned j=0;j<seq.length()-sh.length_+1; ++j) 
// 			{
// 				// std::string window = str_seq.substr(j, w);
// 				// pre_kmer = computeKmerHash(window);
// 				// cout<<"sub seq = "<<window<<endl;
// 				// && prev_key != key
// 				if(sh.set_seed(key, &seq[j]) )//get key via seq
// 				{
// 					it->push(key, seqs.position(i, j), range);
// 					// prev_key = key;
// 				}	
					
// 			}
// 		}
// 		//copy all buffer data into ptr set
// 		it->flush();
// 	}

// 	Ptr_set* build_iterators(const ShapeHistogram &hst) const
// 	{
// 		Ptr_set *iterators = new Ptr_set;
// 		for(unsigned i=0;i<VATConsts::seedp;++i)
// 			(*iterators)[0][i] = ptr_begin(i);

// 		for(unsigned i=1;i<VATConsts::seqp;++i) {
// 			for(unsigned j=0;j<VATConsts::seedp;++j)
// 			{
// 				(*iterators)[i][j] = (*iterators)[i-1][j] + hst[i-1][j];
// 				//cout<<"iterator = "<<(*iterators)[i][j]<<", hst = "<<(hst[i-1][j])<<endl;
// 			}	

// 		}
// 		return iterators;
// 	}
// 	void sortSIMD(Tuple* s ,Tuple* e ) 
// 	// void sortSIMD() 
// 	{
// 		uint64_t N = e - s;
// 		uint64_t p = N % 8;
// 		uint64_t n = N;
// 		std::pair<uint64_t, uint64_t> *rand_arr;
// 		std::pair<uint64_t, uint64_t> *soln_arr;
// 		if(p != 0)
// 		{   
// 			n = N + 8 - p  ;
// 		}
// 		// allocate n memory
// 		aligned_init<std::pair<uint64_t, uint64_t> >(rand_arr, n);
// 		for (size_t i = 0; i < n; i++)
// 		{
// 			if (i < N )
// 			{
// 				uint64_t k = (uint64_t)s[i].key;
// 				uint64_t v = (uint64_t)s[i].value;
// 				// std::pair<uint64_t ,uint64_t> pair_(k,v);
// 				rand_arr[i].first = k;
// 				rand_arr[i].second = v;
// 			}else
// 			{
// 				// std::pair<uint64_t ,uint64_t> Padding(MAX_uint64,MAX_uint64);
// 				rand_arr[i].first = MAX_uint64;
// 				rand_arr[i].second = MAX_uint64;
// 			}	
// 		}
//    		aligned_init<std::pair<uint64_t ,uint64_t> >(soln_arr, n);
// 		std::copy(rand_arr, rand_arr + n, soln_arr);
// 		avx2::SIMDSort(n, soln_arr);
// 		for (size_t i = 0; i < n; i++)
//     	{
// 			if ((soln_arr[i].first != MAX_uint64) && (soln_arr[i].second != MAX_uint64))
// 			{
// 				unsigned k = (unsigned)soln_arr[i].first;
// 				_pos v = (_pos)soln_arr[i].second;
// 				s[i].key = k;
// 				s[i].value = v;
// 			}
//     	}
// 		delete rand_arr;
// 		delete soln_arr;
// 	}

// 	struct Sort_context
// 	{
// 		Sort_context(SortedTuples &sl):
// 			sl (sl)
// 		{ }

// 		bool is_power_of_two(int n) {
// 			return (n != 0) && ((n & (n - 1)) == 0);
// 		}
// 		void operator()(unsigned thread_id ,unsigned seedp) const
// 		{
// 			int n = sl.ptr_end(seedp) - sl.ptr_begin(seedp) ;
// 			if((n >= 8 && ((n != 0) && ((n & (n - 1)) == 0))) && VATParameters::simd_sort)
// 			{
// 				// cout<<"Using SIMD Sort "<<n<<endl;
// 				sl.sortSIMD(sl.ptr_begin(seedp), sl.ptr_end(seedp));
// 			}else
// 			{
				
// 				std::sort(sl.ptr_begin(seedp), sl.ptr_end(seedp)); 
				
// 			}

// 		}
// 		SortedTuples &sl;
// 	};

// 	struct Limits : vector<size_t>
// 	{
// 		Limits(const ShapeHistogram &hst, const seedp_range &range)
// 		{
// 			TimerTools timer ("Computing limits", false);
// 			this->push_back(0);
// 			for(unsigned i=0;i<VATConsts::seedp;++i) {
// #ifdef EXTRA
// 				log_stream << i << ' ' << partition_size(hst, i) << endl;
// #endif
// 				this->push_back(this->operator[](i) + (range.contains(i) ? partition_size(hst, i) : 0));
// 				// cout<<"limit = "<<this->operator[](i)<<endl;
// 			}
// 		}
// 	};

// 	const Limits limits_;
// 	Tuple *data_;

// };

// #endif /* SORTED_LIST_H_ */

