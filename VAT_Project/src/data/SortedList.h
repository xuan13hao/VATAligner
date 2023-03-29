
#ifndef SORTED_LIST_H_
#define SORTED_LIST_H_

#include <memory>
#include <map>
#include "../tools/sort_simd/include/avx256/simd_sort.h"
#include "../tools/util.h"
#include "SeedHistogram.h"
#include "../basic/PackedLocations.h"
#define MAX_uint64 9999999

using std::auto_ptr;
using namespace std;
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
		void set_key(unsigned i)
		{
			this->key = i;
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
		bool operator<(const Tuple &rhs) const
		{ 
			return key < rhs.key; 
		}
		// entry operator++() const
		// { return key < rhs.key; }
		unsigned	key;
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
		launch_scheduled_thread_pool(build_context, VATConsts::seqp, VATParameters::threads());
		timer.go("Sorting seed list");
		Sort_context sort_context (*this);
		launch_scheduled_thread_pool(sort_context, VATConsts::seedp, VATParameters::threads());
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
			const unsigned p (seed_partition(key));
			if(range.contains(p)) {
				assert(n[p] < BUFFER_SIZE);
				// buf[p][n[p]++] = Tuple (key, value);
				buf[p][n[p]++] = Tuple (seed_partition_offset(key), value);
				if(n[p] == BUFFER_SIZE)
					flush(p);
			}
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
		auto_ptr<buffered_iterator> it (new buffered_iterator(ptr));
		for(size_t i=begin;i<end;++i) 
		{
			const sequence<const _val> seq = seqs[i];
			// cout<<"seq = "<<seq<<endl;
			// cout<<"------------ "<<endl;
			if(seq.length()<sh.length_) continue;
	
			for(unsigned j=0;j<seq.length()-sh.length_+1; ++j) 
			{
				// cout<<"seq = "<<seq[j]<<endl;
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
		/*
		int test[64] = {3,1,2,4,3,1,2,4,3,1,2,4,3,1,2,4,3,1,2,4,3,1,2,4,3,1,2,4,3,1,2,4,3,1,2,4,3,1,2,4,3,1,2,4,3,1,2,4,3,1,2,4,3,1,2,4,3,1,2,4,3,1,2,4};
		int* rand_arr;
  		int *soln_arr;
		aligned_init<int>(rand_arr, 64);
		for (size_t i = 0; i < 64; i++)
		{
			rand_arr[i] = test[i];
		}
		// cout<<"============1==========="<<endl;
		aligned_init<int>(soln_arr, 64);
		// cout<<"======================="<<endl;
  		std::copy(rand_arr, rand_arr + 64, soln_arr);
		// cout<<"============2==========="<<endl;
		avx2::SIMDSort(64, soln_arr);
		// cout<<"============3==========="<<endl;
		for (size_t i = 0; i < 64; i++)
		{
			std::cout<<"arr = "<<soln_arr[i]<<std::endl;
		}
		delete soln_arr;
		delete rand_arr;

		*/
		int N = e - s;
		for (size_t i = 0; i < N; i++)
		{
			cout<<"arr "<<i<<"= "<<s[i].key<<":"<<s[i].value<<endl;
			// soln_arr++;
		}
		int p = N % 8;
		int n = N;
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
		cout<<"====================================="<<endl;
    	for (size_t i = 0; i < n; i++)
		{
			cout<<"arr "<<i<<"= "<<rand_arr[i].first<<":"<<rand_arr[i].second<<endl;
		}
   		aligned_init<std::pair<uint64_t ,uint64_t> >(soln_arr, n);
		std::copy(rand_arr, rand_arr + n, soln_arr);
		cout<<"====================================="<<endl;
    	for (size_t i = 0; i < n; i++)
		{
			cout<<"arr "<<i<<"= "<<soln_arr[i].first<<":"<<soln_arr[i].second<<endl;
			// soln_arr++;
		}
		cout<<"===================before=================="<<endl;
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
			// cout<<"=========================1"<<endl;
			// SortedList::const_iterator j = sl.get_partition_cbegin(seedp);
			// while(!j.at_end()) 
			// {
			// 	cout<<"key = "<<j.key()<<endl;
			// 	++j;
			// }
			// cout<<"=========================2"<<endl;
			// Tuple* s1 = sl.ptr_begin(seedp);
			// cout<<"s = "<<s1->key<<endl;
			// Tuple* e1 = sl.ptr_end(seedp);
			// cout<<"e = "<<e1->key<<endl;
			// cout<<"=========================3"<<endl;
			// cout<<"=========================start"<<endl;
			// int c = 0;
			// Tuple* s = sl.ptr_begin(seedp);
			// // cout<<"s = "<<s->key<<endl;
			// Tuple* e = sl.ptr_end(seedp);
			// while (s != e)
			// {
			// 	c = c + 1;
			// 	cout<<"k = "<<s->key<<",v = "<<s->value<<endl;
			// 	s++;
			// }
			// cout<<"c = "<<c<<endl;
			// cout<<"=========================beforesort"<<endl;
			int n = sl.ptr_end(seedp) - sl.ptr_begin(seedp) ;
			// cout<<"n = "<<n<<endl;
			if(n >= 8 && ((n != 0) && ((n & (n - 1)) == 0)))
			{
				sl.sortSIMD(sl.ptr_begin(seedp), sl.ptr_end(seedp));
			}else
			{
				std::sort(sl.ptr_begin(seedp), sl.ptr_end(seedp)); 
				
			}

			// int c1 = 0;
			// std::sort(sl.ptr_begin(seedp), sl.ptr_end(seedp)); 
			// cout<<"=========================aftersort"<<endl;
			// // // int c = 0;
			// Tuple* s1 = sl.ptr_begin(seedp);
			// // cout<<"s = "<<s->key<<endl;
			// Tuple* e1 = sl.ptr_end(seedp);
			// while (s1 != e1)
			// {
			// 	c1 = c1 + 1;
			// 	cout<<"k = "<<s1->key<<",v = "<<s1->value<<endl;
			// 	s1++;
			// }
			// cout<<"c1 = "<<c1<<endl;
			// cout<<"=========================end"<<endl;
			// cout<<"c = "<<c<<endl;
			// cout<<"after sort"<<endl;
			// sl.get_partition_cbegin(seedp);
			// const typename SortedList<_loc>::Type &idx,
			// cout<<"+++++++++++++++++++++++++1"<<endl;
			// SortedList::const_iterator i = sl.get_partition_cbegin(seedp);
			
			// while(!i.at_end()) 
			// {
			// 	cout<<"key = "<<i.key()<<endl;
			// 	++i;
			// }
			// cout<<"+++++++++++++++++++++++++2"<<endl;
			// int c = 0;
			// Tuple* s = sl.ptr_begin(seedp);
			// // cout<<"s = "<<s->key<<endl;
			// Tuple* e = sl.ptr_end(seedp);
			// while (s < e)
			// {
			// 	c = c + 1;
			// 	cout<<"s = "<<s->key<<endl;
			// 	s++;
			// }
			// cout<<"c = "<<c<<endl;
			// cout<<"e = "<<e->key<<endl;
			// cout<<"+++++++++++++++++++++++++3"<<endl;
			// cout<<endl;
			// cout<<"sl.ptr_begin(seedp) = "<<sl.ptr_begin(seedp)<<", sl.ptr_end(seedp) = "<<sl.ptr_end(seedp)<<endl;
			// cout<<"sl.ptr_begin(seedp)= "<<(sl.ptr_begin(seedp))<<endl;
			// cout<<"sl.ptr_begin(seedp)= "<<(sl.ptr_begin(seedp))++<<endl;
			// for(int i = 1; i<VATConsts::seedp;i++)
			// {
			// 	cout<<"sl.ptr_begin(seedp)= "<<(sl.ptr_begin(seedp))[i]<<endl;
			// }
			
			// for(SortedList::iterator it = sl.ptr_begin(seedp); it.at_end();it++)
			// {
			// 	//cout<<"sort list = "<<it<<endl;
			// }
			// while(true)
			// {
			// 	if ((sl.ptr_begin(seedp)++) != sl.ptr_end(seedp))
			// 	{
			// 		cout<<"sl.ptr_begin(seedp)= "<<(sl.ptr_begin(seedp)++)<<endl;	
			// 	}
				
				
			// }
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
