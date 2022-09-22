#ifndef __OBTAINSEEDS_H__
#define __OBTAINSEEDS_H__


#include <vector>
#include <assert.h>
#include <memory>
#include "../utils/merge_sort.h"
#include "../basic/DiagonalSeeds.h"
#include "../search/trace_pt_buffer.h"
#include "../utils/map.h"
#include "../utils/task_queue.h"
#include "../utils/Queue.h"
#include "../basic/Hits.h"
#include "../utils/async_buffer.h"
#include "../basic/Statistics.h"
#include "../search/XdropUngapped.h"
#include "../utils/text_buffer.h"
#include "../out/output_buffer.h"
#include "../utils/binary_file.h"
#include "link_segments.h"
#include "../dp/floating_sw.h"
// #include "./FindSeedsChain.h"
using std::unique_ptr;
// unique_ptr<vector<DiagonalSeeds> >dseed;

using std::vector;
using std::unique_ptr;

struct OutputWriter
{
	OutputWriter(OutputStreamer* f):
		f_ (f)
	{ }

	void operator()(Text_buffer &buf)
	{
		f_->write(buf.get_begin(), buf.size());
		buf.clear();
	}
private:
	OutputStreamer* const f_;
};
static thread_specific_ptr<vector<DiagonalSeeds> > seeds_ptr;
Tls<vector<DiagonalSeeds> > seed_(seeds_ptr);


template<typename _val>
void accessSequence(vector<Segment<_val> > &matches,
		Statistics &stat,
		vector<local_match<_val> > &local,
		unsigned *padding,
		size_t db_letters,
		unsigned dna_len,
		typename Trace_pt_buffer::Vector::iterator &begin,
		typename Trace_pt_buffer::Vector::iterator &end,
		vector<DiagonalSeeds>& diagonalsegment_
		)
{
	string qry,sbj;
	std::sort(begin, end, hit::cmp_normalized_subject);
	const unsigned q_num (begin->query_);
	const sequence<const _val> query (QuerySeqs<_val>::get()[q_num]);
	const unsigned frame = q_num % query_contexts();
	const unsigned query_len = query.length();
	padding[frame] = VATParameters::read_padding<_val>(query_len);
	const SequenceSet<_val> *ref = ReferenceSeqs<_val>::data_;
	for(typename Trace_pt_buffer::Vector::iterator i = begin; i != end; ++i) 
	{
		if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) 
		{
			stat.inc(Statistics::DUPLICATES);
			continue;
		}
		std::pair<size_t, size_t> l = ref->local_position(i->subject_);
		const _val* sbj = ref->data(i->subject_);
		const _val* qry = &query[i->seed_offset_];
		
		// cout<<"i->query_ = "<<(int)i->query_<<", l.second = "<<(int)l.second<<", i->subject = "<<i->subject_<<endl;
		DiagonalSeeds ds = xdrop_ungapped<_val>(qry, sbj,(int)i->seed_offset_,(int)l.second,*i);
		// const _val* sbj1 = ref->data(ds.i+29);
		// const _val* qry1 = &query[ds.j+28];
		// cout<<"subject = "<<AlphabetAttributes<_val>::ALPHABET[*sbj]<<", query = "<<AlphabetAttributes<_val>::ALPHABET[*qry]<<endl;
		cout<<"i = "<<ds.i<<", j = "<<ds.j<<",len= "<<ds.len<<", sbj id ="<<l.first<<", q id ="<<ds.hit_.query_<<endl;
		// seed_->push_back(ds);
		diagonalsegment_.push_back(ds);
	}	
	// for (size_t i = 0; i < diagonalsegment_.size(); i++)
	// {
	// 	cout<<diagonalsegment_[i].qry_id<<"\t"<<diagonalsegment_[i].sbj_id<<"\t"<<diagonalsegment_[i].i<<"\t"<<diagonalsegment_[i].j<<endl;
	// }
}


template<typename _val>
void accessRead(
		Statistics &stat,
		typename Trace_pt_buffer::Vector::iterator &begin,
		typename Trace_pt_buffer::Vector::iterator &end,
		vector<DiagonalSeeds>& ds_
		)
{

	static thread_specific_ptr<vector<local_match<_val> > > local_ptr;
	static thread_specific_ptr<vector<Segment<_val> > > matches_ptr;
	static thread_specific_ptr<vector<char> > transcript_ptr;
	static thread_specific_ptr<vector<DiagonalSeeds> > seeds_ptr;

	Tls<vector<Segment<_val> > > matches (matches_ptr);
	Tls<vector<local_match<_val> > > local (local_ptr);
	Tls<vector<char> > transcript_buf (transcript_ptr);
	Tls<vector<DiagonalSeeds> > seeds_(seeds_ptr);
	local->clear();
	matches->clear();
	transcript_buf->clear();
	seeds_->clear();

	assert(end > begin);
	const size_t hit_count = end - begin;
	local->reserve(hit_count);
	const unsigned contexts = query_contexts();
	const unsigned query = begin->query_/contexts;
	const size_t query_len (QuerySeqs<_val>::data_->length(query*contexts));
	const size_t source_query_len = query_len;
	const size_t db_letters = ref_header.letters;
	unsigned padding[6];
	typedef Map<typename vector<hit >::iterator,typename hit::template Query_id<1> > Map_t;
	Map_t hits_ (begin, end);
	typename Map_t::Iterator i = hits_.begin();
	while(i.valid()) 
	{

		accessSequence<_val>(*matches, stat, *local, padding, db_letters, source_query_len, i.begin(), i.end(),*seeds_);
		++i;
	}
	vector<DiagonalSeeds> ds;
	vector<DiagonalSeeds>::iterator it = seeds_->begin();
	// cout<<"size = "<<seeds_->size()<<endl;
	if (seeds_->size() > 0)
	{
			while (it < seeds_->end())
			{
				
				// ds.push_back(*it);
				// dseed->push_back(*it);
				seed_->push_back(*it);
				// cout<<(*it).qry_id<<"\t"<<(*it).sbj_id<<"\t"<<(*it).i<<"\t"<<(*it).j<<endl;
				++it;
			}
	}
	
	
}

vector<DiagonalSeeds> seeds;
template<typename _val,unsigned _d>
void accessQueries(typename Trace_pt_list::iterator begin,
		typename Trace_pt_list::iterator end,
		Statistics &st,
		vector<DiagonalSeeds>& ds_
		)
{
	typedef Map<typename vector<hit>::iterator,typename hit::template Query_id<_d> > Map_t;
	Map_t hits_ (begin, end);
	typename Map_t::Iterator i = hits_.begin();
	while(i.valid()&& !exception_state()) 
	{
		accessRead<_val>(st, i.begin(), i.end(),ds_);

		++i;

	}
	// cout<<"size = "<<ds_.size()<<endl;
}

template<typename _val, typename _buffer>
struct AlignContext
{
	AlignContext(Trace_pt_list &trace_pts, OutputStreamer* output_file):
		trace_pts (trace_pts),
		output_file (output_file),
		writer (output_file),
		queue (VATParameters::threads()*8, writer)
	{ }
	void operator()(unsigned thread_id)
	{
		
		Statistics st;
		size_t i=0;
		typename Trace_pt_list::Query_range query_range (trace_pts.get_range());
		_buffer *buffer = 0;
		// seeds_->clear();
		while(queue.get(i, buffer, query_range) && !exception_state()) 
		{
			try 
			{

				accessQueries<_val,1>(query_range.begin, query_range.end, st, ds1);
				// cout<<"==========================================================================="<<endl;
				queue.push(i);
			} catch(std::exception &e) 
			{
				exception_state.set(e);
				queue.wake_all();
			}
			
		}
		statistics += st;
		// for (size_t i = 0; i < ds1.size(); i++)
		// {
		// 		cout<<ds1[i].qry_id<<"\t"<<ds1[i].sbj_id<<"\t"<<ds1[i].i<<"\t"<<ds1[i].j<<"\t"<<ds1[i].len<<endl;
		// }
	}
	vector<DiagonalSeeds> ds1;
	Trace_pt_list &trace_pts;
	OutputStreamer* output_file;
	OutputWriter writer;
	Task_queue<_buffer,OutputWriter> queue;
	
};

template<typename _val, typename _locr, typename _locl>
void accessQueries(Trace_pt_buffer &trace_pts, OutputStreamer* output_file)
{
	Trace_pt_list v;
	for(unsigned bin=0;bin<trace_pts.bins();++bin) 
	{
		TimerTools timer ("Loading trace points", false);
		trace_pts.load(v, bin);
		merge_sort(v.begin(), v.end(), VATParameters::threads());
		v.init();
		timer.go("Computing Seeds");
		if(ref_header.n_blocks > 1) 
		{
			// accessQueries<_val,1>(query_range.begin, query_range.end, *buffer, st, ds);
			AlignContext<_val,Temp_output_buffer<_val> > context (v, output_file);
			launch_thread_pool(context, VATParameters::threads());
			vector<DiagonalSeeds> ds = context.ds1;
			for (size_t i = 0; i < ds.size(); i++)
			{
				cout<<ds[i].qry_id<<"\t"<<ds[i].sbj_id<<"\t"<<ds[i].i<<"\t"<<ds[i].j<<"\t"<<ds[i].len<<endl;
			}
		} else 
		{
			// accessQueries<_val,1>(query_range.begin, query_range.end, *buffer, st, ds);
			AlignContext<_val,Output_buffer<_val> > context (v, output_file);
			launch_thread_pool(context, VATParameters::threads());
			// cout<<"==========================================================================="<<endl;
			// vector<DiagonalSeeds> ds = context.ds1;
			// for (size_t i = 0; i < ds.size(); i++)
			// {
			// 	cout<<ds[i].qry_id<<"\t"<<ds[i].sbj_id<<"\t"<<ds[i].i<<"\t"<<ds[i].j<<"\t"<<ds[i].len<<endl;
			// }
		}
		
	}
	// vector<DiagonalSeeds>::iterator i = seed_->begin();
	// cout<<"size = "<<seed_->size()<<endl;
	// while (i<seed_->end())
	// {
	// 	cout<<i->qry_id<<"\t"<<i->sbj_id<<endl;
	// 	i++;
	// }
	
	

		// for (size_t i = 0; i < ds_total.size(); i++)
		// {
		// 		cout<<ds_total[i]->qry_id<<"\t"<<ds_total[i].sbj_id<<"\t"<<ds_total[i].i<<"\t"<<ds_total[i].j<<"\t"<<ds_total[i].len<<endl;
		// }
		// cout<<"size = "<<seed_->size()<<endl;
}
/*
0	0	1	25	24
1	0	1	25	24
2	0	1	25	24

0	0	1	25
1	0	1	25
2	0	1	25
0	1	0	16
1	1	0	16
2	1	0	16
0	2	0	33
1	2	0	33
2	2	0	33


*/

#endif // __OBTAINSEEDS_H__