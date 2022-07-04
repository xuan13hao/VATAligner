

#ifndef TRACE_PT_BUFFER_H_
#define TRACE_PT_BUFFER_H_

#include "../tools/async_buffer.h"
#include "../basic/Hits.h"

using std::auto_ptr;

template<typename _locr, typename _locl>
struct Trace_pt_buffer : public AsynchronousBuffer<Hits<_locr,_locl> >
{
	Trace_pt_buffer(size_t input_size, const string &tmpdir, bool mem_buffered):
		AsynchronousBuffer<Hits<_locr,_locl> > (input_size, tmpdir, mem_buffered ? mem_bins : file_bins)
	{ 
		
	}
	enum { mem_bins = 1, file_bins = 4 };
	static Trace_pt_buffer *instance;
};

template<typename _locr, typename _locl> Trace_pt_buffer<_locr,_locl>* Trace_pt_buffer<_locr,_locl>::instance;

template<typename _locr, typename _locl>
class Trace_pt_list : public vector<Hits<_locr,_locl> >
{
	public:
	void init()
	{
		pos_ = this->begin();
		total_ = 0;
		count_ = 1;
#ifdef PRE_PARTITION
		p_.clear();
		p_.push_back(0);
		idx_ = 0;
		const unsigned c = query_contexts();
		typename vector<hit<_locr,_locl> >::iterator i = this->begin();
		unsigned total=0,count=1;
		for(; i < this->end();) {
			unsigned n=0;
			const unsigned min_size = std::max(4*total/count/5 + 1, program_options::fetch_size);
			for(;i<this->end() && n<min_size;) {
				const unsigned q = i->query_/c;
				for(; i<this->end() && i->query_/c == q; ++i)
					++n;
			}
			++count;
			total += n;
			p_.push_back(i - this->begin());
		}
		p_.push_back(i - this->begin());
#endif
	}
	struct Query_range
	{
		Query_range(Trace_pt_list &parent):
			parent_ (parent)
		{ 
			
		}
#ifndef PRE_PARTITION
		bool operator()()
		{

			begin = parent_.pos_;
			//end = std::min(std::max(begin + 3*parent_.total_/parent_.count_/4 + 1, begin+program_options::fetch_size), parent_.end());
			end = std::min(begin + 3*parent_.total_/parent_.count_/4 + 1, parent_.end());
			if(end >= parent_.end())
				return false;
			const unsigned c = query_contexts(), q = end->query_/c;
			for(; end<parent_.end() && end->query_/c == q; ++end);
			parent_.pos_ = end;
			parent_.total_ += end - begin;
			++parent_.count_;
			return end < parent_.end();
		}
#else
		bool operator()()
		{
			begin = parent_.begin()+parent_.p_[parent_.idx_];
			end = parent_.begin()+parent_.p_[parent_.idx_+1];
			printf("%lu %lu %lu\n", parent_.p_[parent_.idx_], parent_.p_[parent_.idx_+1], parent_.p_[parent_.idx_+1]-parent_.p_[parent_.idx_]);
			++parent_.idx_;
			return parent_.idx_ < parent_.p_.size()-1;
		}
#endif
		typename Trace_pt_list::iterator begin, end;
	private:
		Trace_pt_list &parent_;
	};
	Query_range get_range()
	{ return Query_range (*this); }
private:
	typename vector<Hits<_locr,_locl> >::iterator pos_;
#ifdef PRE_PARTITION
	vector<size_t> p_;
	unsigned idx_;
#else
	size_t total_, count_;
#endif
};

#endif /* TRACE_PT_BUFFER_H_ */

