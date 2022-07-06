
#ifndef SORTED_LIST_H_
#define SORTED_LIST_H_

#include "../tools/util.h"
#include "SeedHistogram.h"
#include "../basic/PackedLocations.h"
#include <memory>
using std::auto_ptr;

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
		bool operator<(const Tuple &rhs) const
		{ 
			return key < rhs.key; 
		}
		// entry operator++() const
		// { return key < rhs.key; }
		unsigned	key;
		unsigned 	seed2;
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

	struct Sort_context
	{
		Sort_context(SortedList &sl):
			sl (sl)
		{ }
		void operator()(unsigned thread_id ,unsigned seedp) const
		{
			// SortedList::const_iterator j = sl.get_partition_cbegin(seedp);
			// while(!j.at_end()) 
			// {
			// 	cout<<"key = "<<j.key()<<endl;
			// 	++j;
			// }
			
			std::sort(sl.ptr_begin(seedp), sl.ptr_end(seedp)); 
			// cout<<"after sort"<<endl;
			// sl.get_partition_cbegin(seedp);
			// const typename SortedList<_loc>::Type &idx,
			// SortedList::const_iterator i = sl.get_partition_cbegin(seedp);
			// while(!i.at_end()) 
			// {
			// 	cout<<"key = "<<i.key()<<endl;
			// 	++i;
			// }
			// cout<<endl;
			// entry* s = sl.ptr_begin(seedp);
			// cout<<"s = "<<s->key<<endl;
			// entry* e = sl.ptr_end(seedp);
			// cout<<"e = "<<e->key<<endl;
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
