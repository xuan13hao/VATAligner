

#ifndef ALIGN_QUERIES_H_
#define ALIGN_QUERIES_H_

#include "../util/merge_sort.h"
#include "../search/trace_pt_buffer.h"
#include "../util/map.h"
#include "align_read.h"
#include "../util/task_queue.h"

using std::vector;

struct Output_writer
{
	Output_writer(Output_stream* f):
		f_ (f)
	{ }
	void operator()(Text_buffer &buf)
	{
		f_->write(buf.get_begin(), buf.size());
		buf.clear();
	}
private:
	Output_stream* const f_;
};

template<typename _val, typename _locr, typename _locl, unsigned _d>
void align_queries(typename Trace_pt_list<_locr,_locl>::iterator begin,
		typename Trace_pt_list<_locr,_locl>::iterator end,
		Output_buffer<_val> &buffer,
		Statistics &st)
{
	typedef Map<typename vector<hit<_locr,_locl> >::iterator,typename hit<_locr,_locl>::template Query_id<_d> > Map_t;
	Map_t hits (begin, end);
	typename Map_t::Iterator i = hits.begin();
	while(i.valid() && !exception_state()) {
		cout<<"align_read 1"<<endl;
		align_read<_val,_locr,_locl>(buffer, st, i.begin(), i.end());
		++i;
	}
}

template<typename _val, typename _locr, typename _locl, typename _buffer>
struct Align_context
{
	Align_context(Trace_pt_list<_locr,_locl> &trace_pts, Output_stream* output_file):
		trace_pts (trace_pts),
		output_file (output_file),
		writer (output_file),
		queue (VATParameters::threads()*8, writer)
	{ }
	void operator()(unsigned thread_id)
	{
		Statistics st;
		size_t i=0;
		typename Trace_pt_list<_locr,_locl>::Query_range query_range (trace_pts.get_range());
		_buffer *buffer = 0;
		while(queue.get(i, buffer, query_range) && !exception_state()) {
			try {

				cout<<"Align_context"<<endl;
				//switch(1) {
				// case 6:
				// 	align_queries<_val,_locr,_locl,6>(query_range.begin, query_range.end, *buffer, st);
				// 	break;
				// case 2:
				// 	align_queries<_val,_locr,_locl,2>(query_range.begin, query_range.end, *buffer, st);
				// 	break;
				//default:
					align_queries<_val,_locr,_locl,1>(query_range.begin, query_range.end, *buffer, st);
									cout<<"Align_context 1"<<endl;

			//	}
				queue.push(i);
			} catch(std::exception &e) {
				exception_state.set(e);
				queue.wake_all();
			}
		}
		statistics += st;
		// cout<<"Align_context"<<endl;
	}
	Trace_pt_list<_locr,_locl> &trace_pts;
	Output_stream* output_file;
	Output_writer writer;
	Task_queue<_buffer,Output_writer> queue;
};

template<typename _val, typename _locr, typename _locl>
void align_queries(const Trace_pt_buffer<_locr,_locl> &trace_pts, Output_stream* output_file)
{
	Trace_pt_list<_locr,_locl> v;
	for(unsigned bin=0;bin<trace_pts.bins();++bin) 
	{
		log_stream << "Processing query bin " << bin+1 << '/' << trace_pts.bins() << '\n';
		task_timer timer ("Loading trace points", false);
		trace_pts.load(v, bin);
		timer.go("Sorting trace points");
		merge_sort(v.begin(), v.end(), VATParameters::threads());
		v.init();
		timer.go("Computing alignments");
		if(ref_header.n_blocks > 1) {
			Align_context<_val,_locr,_locl,Temp_output_buffer<_val> > context (v, output_file);
			launch_thread_pool(context, VATParameters::threads());
		} else {
			Align_context<_val,_locr,_locl,Output_buffer<_val> > context (v, output_file);
			launch_thread_pool(context, VATParameters::threads());
		}
	}
	cout<<"align_queries"<<endl;
}

#endif /* ALIGN_QUERIES_H_ */
