
#ifndef FREQUENCY_MASKING_H_
#define FREQUENCY_MASKING_H_

#include <vector>

using std::vector;
using std::auto_ptr;
// using std::unique_ptr;
template<typename _val>
class Masked_sequence_set : public SequenceSet<_val>
{
	public:

	Masked_sequence_set(Input_stream &file):
			SequenceSet<_val> (file)
	{ }

	template<typename _loc>
	void build_masking(unsigned sid, const seedp_range &range, typename SortedList<_loc>::Type &idx)
	{
		TimerTools timer ("Counting low complexity seeds", false);
		vector<unsigned> counts (VATConsts::seedp);
		Count_context<_loc> count_context (idx, counts);
		launch_scheduled_thread_pool(count_context, VATConsts::seedp, VATParameters::threads());

		timer.finish();
		size_t n = 0;
		for(unsigned i=range.begin();i<range.end();++i) {
			n += counts[i];
			const size_t ht_size (std::max(static_cast<size_t>(static_cast<float>(counts[i]) * 1.3), static_cast<size_t>(counts[i] + 1)));
			pos_filters[sid][i] = auto_ptr<filter_table> (new filter_table(ht_size));
		}
		log_stream << "Hit cap = " << VATParameters::hit_cap << std::endl;
		log_stream << "Low complexity seeds = " << n << std::endl;

		timer.go("Building position filter");
		Build_context<_loc> build_context(idx, sid, counts, *this);
		launch_scheduled_thread_pool(build_context, VATConsts::seedp, VATParameters::threads());
		timer.finish();
		log_stream << "Masked positions = " << std::accumulate(counts.begin(), counts.end(), 0) << std::endl;
	}

	bool get_masking(const _val *pos, unsigned sid) const
	{
		uint64_t seed;
		// cout<<"get_masking 1"<<endl;
		ShapeConfigures::get().get_shape(sid).set_seed(seed, pos);
		// cout<<"get_masking 2"<<endl;
		const filter_table::Tuple *e;
		// cout<<"get_masking 3"<<endl;
		// size_t s = pos_filters[sid][seed_partition(seed)]->size();
		// cout<<"get_masking 5"<<endl;
		if((e = pos_filters[sid][seed_partition(seed)]->operator [](seed_partition_offset(seed))) != 0) 
		{
					
			const size_t offset (pos - this->data(0));
					// cout<<"get_masking 4"<<endl;
			return !position_filter(offset, e->value, seed_partition_offset(seed));
		} else
			return false;
	}

	virtual ~Masked_sequence_set()
	{ }

private:

	template<typename _loc>
	struct Count_context
	{
		Count_context(const typename SortedList<_loc>::Type &idx, vector<unsigned> &counts):
			idx (idx),
			counts (counts)
		{ }
		void operator()(unsigned thread_id, unsigned seedp) const
		{
			unsigned n = 0;
			typename SortedList<_loc>::Type::const_iterator i = idx.get_partition_cbegin(seedp);
			while(!i.at_end()) {
				if(i.n > VATParameters::hit_cap)
					++n;
				++i;
			}
			counts[seedp] = n;
		}
		const typename SortedList<_loc>::Type &idx;
		vector<unsigned> &counts;
	};

	template<typename _loc>
	struct Build_context
	{
		Build_context(const typename SortedList<_loc>::Type &idx, unsigned sid, vector<unsigned> &counts, Masked_sequence_set<_val> &seqs):
			idx (idx),
			sid (sid),
			counts (counts),
			seqs (seqs)
		{ }
		void operator()(unsigned thread_id, unsigned seedp)
		{
			unsigned n = 0;
			typename SortedList<_loc>::Type::iterator i = idx.get_partition_begin(seedp);
			while(!i.at_end()) {
				if(i.n > VATParameters::hit_cap)
					n += seqs.mask_seed_pos<_loc>(i, sid, seedp);
				++i;
			}
			counts[seedp] = n;
		}
		const typename SortedList<_loc>::Type &idx;
		const unsigned sid;
		vector<unsigned> &counts;
		Masked_sequence_set<_val> &seqs;
	};

	template<typename _loc>
	unsigned mask_seed_pos(typename SortedList<_loc>::Type::iterator &i, unsigned sid, unsigned p)
	{
		const unsigned treshold (filter_treshold(i.n));
		unsigned count (0), k (0);
		for(unsigned j=0;j<i.n;++j)
			if(!position_filter(i[j], treshold, i.key())) {
				mask_seed_pos<_loc>(i[j]);
				++count;
			} else
				*(i.get(k++)) = *(i.get(j));
		i.get(k)->value = 0;
		pos_filters[sid][p]->insert(i.key(), treshold);
		return count;
	}

	template<typename _loc>
	void mask_seed_pos(_loc pos)
	{
		_val *x = this->data(pos);
		*x = set_critical(*x);
	}

private:

	typedef hash_table<uint32_t, uint8_t, value_compare<uint8_t, 0>, murmur_hash> filter_table;
	auto_ptr<filter_table> pos_filters[VATConsts::max_shapes][VATConsts::seedp];

};

#endif /* FREQUENCY_MASKING_H_ */
