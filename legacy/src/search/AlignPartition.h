

#ifndef ALIGN_RANGE_H_
#define ALIGN_RANGE_H_

#include "Align.h"
#include "../basic/Statistics.h"

template<typename _val, typename _locr, typename _locq, typename _locl>
void align_range(_locq q_pos,
				 const typename SortedList<_locr>::Type::const_iterator &s,
				 Statistics &stats,
				 typename Trace_pt_buffer<_locr,_locl>::Iterator &out,
				 unsigned sid)
{
	unsigned i = 0, n=0;

	const _val* query = QuerySeqs<_val>::data_->data(q_pos);
	HitFilter<_val,_locr,_locq,_locl> hf (stats, q_pos, out);
	if(s.n <= VATParameters::hit_cap) 
	{
		stats.inc(Statistics::SEED_HITS, s.n);
		while(i < s.n) 
		{
			align<_val,_locr,_locq,_locl>(q_pos, query, s[i], stats, sid, hf);
			++i;
		}
	} 
	else 
	{
		while(i < s.n && s[i] != 0) 
		{
			assert(position_filter(s[i], filter_treshold(s.n), s.key()));
			align<_val,_locr,_locq,_locl>(q_pos, query, s[i], stats, sid, hf);
			stats.inc(Statistics::SEED_HITS);
			++i;
			++n;
		}
	}
#ifdef EXTRA
	//if(n > 64)
		//printf("%u\n",n);
#endif

	hf.finish();
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void align_range(const typename SortedList<_locq>::Type::const_iterator &q,
				 const typename SortedList<_locr>::Type::const_iterator &s,
				 Statistics &stats,
				 typename Trace_pt_buffer<_locr,_locl>::Iterator &out,
				 const unsigned sid)
{
#ifdef EXTRA
	//if(q.n > 4096)
		//printf("%lu %lu\n",q.n,s.n);
#endif

//maximum number of hits to consider for one seed
	for(unsigned i=0;i<q.n; ++i)
		align_range<_val,_locr,_locq,_locl>(_locq(q[i]), s, stats, out, sid);
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void alignPartition(unsigned hp,
		Statistics &stats,
		unsigned sid,
		typename SortedList<_locr>::Type::const_iterator i,
		typename SortedList<_locq>::Type::const_iterator j,
		unsigned thread_id)
{
	typename Trace_pt_buffer<_locr,_locl>::Iterator* out = new typename Trace_pt_buffer<_locr,_locl>::Iterator (*Trace_pt_buffer<_locr,_locl>::instance, thread_id);
	while(!i.at_end() && !j.at_end() && !exception_state()) 
	{
		if(i.key() < j.key()) 
		{
			++i;
		} else if(j.key() < i.key()) 
		{
			++j;
		} 
		else 
		{
			align_range<_val,_locr,_locq,_locl>(j, i, stats, *out, sid);
			++i;
			++j;
		}
	}
	delete out;
}

#endif /* ALIGN_RANGE_H_ */
