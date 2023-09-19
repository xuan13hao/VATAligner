
#ifndef VIEW_H_
#define VIEW_H_

const unsigned view_buf_size = 32;

class View_writer
{
	public:
	View_writer():
		f_ (VATParameters::output_file + (VATParameters::compression==1?".gz":""),
				VATParameters::compression==1,
				VATParameters::output_file.length()==0 ? OutputStreamer::stdout_sink : OutputStreamer::file_sink)
	{ }
	void operator()(Text_buffer &buf)
	{
		f_.write(buf.get_begin(), buf.size());
		buf.clear();
	}
	OutputStreamer f_;
};

struct View_fetcher
{
	View_fetcher(VATFile &daa):
		daa (daa)
	{ }
	bool operator()()
	{
		n = 0;
		for(unsigned i=0;i<view_buf_size;++i)
			if(!daa.read_query_buffer(buf[i]))
				return false;
			else
				++n;
		return true;
	}
	Binary_buffer buf[view_buf_size];
	unsigned n;
	VATFile &daa;
};

template<typename _val>
class View_context
{
	public:
	View_context(VATFile &daa, View_writer &writer, const Output_format<_val> &format):
		daa (daa),
		writer (writer),
		queue (3*VATParameters::threads(), writer),
		format (format)
	{ 
		
	}
	void operator()(unsigned thread_id)
	{
		try {
			size_t n;
			View_fetcher query_buf (daa);
			Text_buffer *buffer = 0;
			
			while(!exception_state() && queue.get(n, buffer, query_buf)) {

				for(unsigned j=0;j<query_buf.n;++j) {
					VATQueryRecord<_val> r (daa, query_buf.buf[j]);
					// cout<<"format print 1"<<endl;
					for(typename VATQueryRecord<_val>::Match_iterator i = r.begin(); i.good(); ++i) 
					{
						// cout<<"format print 2" <<endl;

						if(i->frame > 2 && VATParameters::forwardonly)
							continue;
						// cout<<"format print 3"<<endl;

						format.print_match(*i, *buffer);

						// cout<<"format print 4"<<endl;
					}
				}
				queue.push(n);
			}
		} catch(std::exception &e) {
			exception_state.set(e);
			queue.wake_all();
		}
	}
	VATFile &daa;
	View_writer &writer;
	Task_queue<Text_buffer,View_writer> queue;
	const Output_format<_val> &format;
};

template<typename _val>
void view(VATFile &daa)
{

	ScoreMatrix::instance = auto_ptr<ScoreMatrix> (new ScoreMatrix(daa.score_matrix(),
					daa.gap_open_penalty(),
					daa.gap_extension_penalty(),
					daa.match_reward(),
					daa.mismatch_penalty(),
					_val ()));

	// cout << "Build version = " << daa.diamond_build() << endl;
	// cout << "DB sequences = " << daa.db_seqs() << endl;
	// cout << "DB sequences used = " << daa.db_seqs_used() << endl;
	// cout << "DB letters = " << daa.db_letters() << endl;

	View_writer writer;
	const Output_format<_val>& format (get_output_format<_val>());
	format.print_header(writer.f_);

	View_context<_val> context (daa, writer, format);
	// cout<<"View_context..."<<endl;

	launch_thread_pool(context, VATParameters::threads());
// cout<<"View_context..."<<endl;

}

void view()
{
	VATFile daa (VATParameters::daa_file);
	if(daa.mode() == dna)
	{	
		VATParameters::algn_type = VATParameters::dna; 
		view<DNA>(daa);
	}
	else
	{
		VATParameters::algn_type = VATParameters::protein; 
		view<Protein>(daa);

	}
		
}

#endif /* VIEW_H_ */
