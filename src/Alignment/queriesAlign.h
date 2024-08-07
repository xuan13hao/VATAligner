

#ifndef ALIGN_QUERIES_H_
#define ALIGN_QUERIES_H_

#include "../tools/merge_sort.h"
#include "../Filter/trace_pt_buffer.h"
#include "../tools/map.h"
#include "readsAlign.h"
#include "../tools/task_queue.h"

using std::vector;

struct Output_writer
{
	Output_writer(OutputStreamer* f):
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

template<typename _val, typename _locr, typename _locl, unsigned _d>
void alignQueries(typename Trace_pt_list<_locr,_locl>::iterator begin,
		typename Trace_pt_list<_locr,_locl>::iterator end,
		Output_buffer<_val> &buffer,
		Statistics &st)
{
	typedef Map<typename vector<Hits<_locr,_locl> >::iterator,typename Hits<_locr,_locl>::template Query_id<_d> > Map_t;
	Map_t hits (begin, end);
	typename Map_t::Iterator i = hits.begin();
	while(i.valid() && !exception_state()) 
	{
		// cout<<"align_read ............"<<endl;
		align_read<_val,_locr,_locl>(buffer, st, i.begin(), i.end());
		++i;
	}
}

template<typename _val, typename _locr, typename _locl, typename _buffer>
struct Align_context
{
	Align_context(Trace_pt_list<_locr,_locl> &trace_pts, OutputStreamer* output_file):
		trace_pts (trace_pts),
		output_file (output_file),
		writer (output_file),
		queue (VATParameters::threads()*8, writer)
	{ }
	void operator()(unsigned thread_id)
	{
		Statistics st;
		size_t i=0;
		//vector to store
		// cout<<"Align_context............."<<endl;
		typename Trace_pt_list<_locr,_locl>::Query_range query_range (trace_pts.get_range());
		_buffer *buffer = 0;
		while(queue.get(i, buffer, query_range) && !exception_state()) {
			try {
				switch(query_contexts()) {
				case 6:
					alignQueries<_val,_locr,_locl,6>(query_range.begin, query_range.end, *buffer, st);
					break;
				default:
					alignQueries<_val,_locr,_locl,1>(query_range.begin, query_range.end, *buffer, st);
				}
					// alignQueries<_val,_locr,_locl,1>(query_range.begin, query_range.end, *buffer, st);

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
	OutputStreamer* output_file;
	Output_writer writer;
	Task_queue<_buffer,Output_writer> queue;
};

template<typename _val, typename _locr, typename _locl>
void alignQueries(const Trace_pt_buffer<_locr,_locl> &trace_pts, OutputStreamer* output_file)
{
	Trace_pt_list<_locr,_locl> v;
	for(unsigned bin=0;bin<trace_pts.bins();++bin) 
	{
		// cout << "Processing query bin " << bin+1 << '/' << trace_pts.bins() << '\n';
		TimerTools timer ("Loading trace points", false);
		trace_pts.load(v, bin);

		merge_sort(v.begin(), v.end(), VATParameters::threads());

		
		v.init();
		timer.go("Generating seeds...");
		if(ref_header.n_blocks > 1) {
			Align_context<_val,_locr,_locl,Temp_output_buffer<_val> > context (v, output_file);
			launch_thread_pool(context, VATParameters::threads());
		} else {
			Align_context<_val,_locr,_locl,Output_buffer<_val> > context (v, output_file);
			launch_thread_pool(context, VATParameters::threads());
		}
	}
	// cout<<"align_queries"<<endl;
}

#endif /* ALIGN_QUERIES_H_ */
