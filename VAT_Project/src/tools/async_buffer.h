
#ifndef ASYNC_BUFFER_H_
#define ASYNC_BUFFER_H_

#include <vector>
#include <boost/ptr_container/ptr_vector.hpp>

namespace io = boost::iostreams;

using std::vector;
using std::string;
using std::endl;
using boost::ptr_vector;

struct Buffer_file_read_exception : public VATException
{
	Buffer_file_read_exception(const char* file_name, size_t count, size_t n):
		VATException (string("Error reading buffer file ") + file_name + " (" + boost::lexical_cast<string>(count) + '/' + boost::lexical_cast<string>(n) + ')')
	{ 
		
	}
};

const unsigned async_buffer_max_bins = 4;

template<typename _t>
class AsynchronousBuffer
{
	public:
	typedef vector<_t> Vector;

	AsynchronousBuffer(size_t input_count, const string &tmpdir, unsigned bins):
		bins_ (bins),
		bin_size_ ((input_count + bins_ - 1) / bins_)
	{
		log_stream << "Async_buffer() " << input_count << ',' << bin_size_ << endl;
		for(unsigned j=0;j<VATParameters::threads();++j)
			for(unsigned i=0;i<bins;++i) {
				tmp_file_.push_back(TempFile ());
				out_.push_back(new OutputStreamer (tmp_file_.back()));
				size_.push_back(0);
			}
	}

	struct Iterator
	{
		Iterator(AsynchronousBuffer &parent, unsigned thread_num):
			parent_ (parent),
			thread_num_ (thread_num)
		{
			for(unsigned i=0;i<parent.bins_;++i) {
				size_[i] = 0;
				out_[i] = parent.get_out(thread_num_, i);
			}
		}
		void push(const _t &x)
		{
			const unsigned bin = x / parent_.bin_size_;
			assert(bin < parent_.bins());
			buffer_[bin*buffer_size+(size_[bin]++)] = x;
			if(size_[bin] == buffer_size)
				flush(bin);
		}
		void flush(unsigned bin)
		{
			out_[bin]->write(&buffer_[bin*buffer_size], size_[bin]);
			parent_.add_size(thread_num_, bin, size_[bin]);
			size_[bin] = 0;
		}
		~Iterator()
		{
			for(unsigned bin=0;bin<parent_.bins_;++bin)
				flush(bin);
		}
	private:
		enum { buffer_size = 65536 };
		_t buffer_[async_buffer_max_bins*buffer_size];
		size_t size_[async_buffer_max_bins];
		OutputStreamer* out_[async_buffer_max_bins];
		AsynchronousBuffer &parent_;
		const unsigned thread_num_;
	};

	void close()
	{
		for(ptr_vector<OutputStreamer>::iterator i=out_.begin();i!=out_.end();++i)
			i->close();
		out_.clear();
		log_stream << "Async_buffer.close() " << endl;
	}

	void load(vector<_t> &data, unsigned bin) const
	{
		size_t size = 0;
		for(unsigned i=0;i<VATParameters::threads();++i)
			size += size_[i*bins_+bin];
		log_stream << "Async_buffer.load() " << size << "(" << (double)size*sizeof(_t)/(1<<30) << " GB)" << endl;
		data.resize(size);
		_t* ptr = data.data();
		for(unsigned i=0;i<VATParameters::threads();++i) {
			Input_stream f (tmp_file_[i*bins_+bin]);
			const size_t s = size_[i*bins_+bin];
			const size_t n = f.read(ptr, s);
			ptr += s;
			f.close();
			if(n != s)
				throw Buffer_file_read_exception(f.file_name.c_str(), s, n);
		}
	}

	unsigned bins() const
	{ 
		return bins_; 
	}

private:

	OutputStreamer* get_out(unsigned threadid, unsigned bin)
	{ 
		return &out_[threadid*bins_+bin]; 
	}

	void add_size(unsigned thread_id, unsigned bin, size_t n)
	{ 
		size_[thread_id*bins_+bin] += n; 
	}

	const unsigned bins_, bin_size_;
	ptr_vector<OutputStreamer> out_;
	vector<size_t> size_;
	vector<TempFile> tmp_file_;

};

#endif /* ASYNC_BUFFER_H_ */
