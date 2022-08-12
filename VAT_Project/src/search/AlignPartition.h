

#ifndef ALIGN_RANGE_H_
#define ALIGN_RANGE_H_

#include "align.h"
#include "../basic/Statistics.h"

template<typename _val, typename _locr, typename _locq, typename _locl>
void alignQueryRange(const typename SortedList<_locq>::Type::const_iterator &q,
				 const typename SortedList<_locr>::Type::const_iterator &s,
				 Statistics &stats,
				 typename Trace_pt_buffer::Iterator &out,
				 const unsigned sid)
{

// printf("%lu %lu\n",q.n,s.n);
//maximum number of hits to consider for one seed
	for(unsigned i=0;i<q.n; ++i)
	{
		_locq q_pos = _locq(q[i]);
		unsigned num = 0, n=0;
		const _val* query = QuerySeqs<_val>::data_->data(q_pos);
		HitFilter<_val,_locr,_locq,_locl> hf (stats, q_pos, out);
		if(s.n <= VATParameters::hit_cap) 
		{
			stats.inc(Statistics::SEED_HITS, s.n);
			while(num < s.n) 
			{
				align<_val,_locr,_locq,_locl>(q_pos, query, s[num], stats, sid, hf);
				++num;
			}
		} 
		else 
		{
			while(num < s.n && s[num] != 0) 
			{
					assert(position_filter(s[num], filter_treshold(s.n), s.key()));
					align<_val,_locr,_locq,_locl>(q_pos, query, s[num], stats, sid, hf);
					stats.inc(Statistics::SEED_HITS);
					++num;
					++n;
			}
		}

		hf.finish();

	}
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void alignPartition(unsigned hp,
		Statistics &stats,
		unsigned sid,
		typename SortedList<_locr>::Type::const_iterator i,
		typename SortedList<_locq>::Type::const_iterator j,
		unsigned thread_id)
{
	typename Trace_pt_buffer::Iterator* out = new typename Trace_pt_buffer::Iterator (*Trace_pt_buffer::instance, thread_id);
	while(!i.at_end() && !j.at_end() && !exception_state()) 
	{
		// cout<<"i key = "<<i.key()<<", j key = "<<j.key()<<endl;
		if(i.key() < j.key()) 
		{
			++i;
		} else if(j.key() < i.key()) 
		{
			++j;
		} 
		else 
		{
			alignQueryRange<_val,_locr,_locq,_locl>(j, i, stats, *out, sid);
			++i;
			++j;
		}
	}
	delete out;
}

#endif /* ALIGN_RANGE_H_ */
