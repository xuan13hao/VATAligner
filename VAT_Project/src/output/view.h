
#ifndef VIEW_H_
#define VIEW_H_

const unsigned view_buf_size = 32;

struct View_writer
{
	View_writer():
		f_ (VATParameters::output_file + (VATParameters::compression==1?".gz":""),
				VATParameters::compression==1,
				VATParameters::output_file.length()==0 ? Output_stream::stdout_sink : Output_stream::file_sink)
	{ }
	void operator()(Text_buffer &buf)
	{
		f_.write(buf.get_begin(), buf.size());
		buf.clear();
	}
	Output_stream f_;
};

struct View_fetcher
{
	View_fetcher(DAA_file &daa):
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
	DAA_file &daa;
};

template<typename _val>
struct View_context
{
	View_context(DAA_file &daa, View_writer &writer, const Output_format<_val> &format):
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
					DAA_query_record<_val> r (daa, query_buf.buf[j]);
					for(typename DAA_query_record<_val>::Match_iterator i = r.begin(); i.good(); ++i) {
						if(i->frame > 2 && VATParameters::forwardonly)
							continue;
						format.print_match(*i, *buffer);
						

					}
				}
				queue.push(n);
			}
		} catch(std::exception &e) {
			exception_state.set(e);
			queue.wake_all();
		}
	}
	DAA_file &daa;
	View_writer &writer;
	Task_queue<Text_buffer,View_writer> queue;
	const Output_format<_val> &format;
};

template<typename _val>
void view(DAA_file &daa)
{

	// cout<<""<<daa.score_matrix()<<" " <<daa.gap_open_penalty()<<" "<<daa.gap_extension_penalty()<<" "<<daa.match_reward()<<" "<<daa.mismatch_penalty()<<endl;
	score_matrix::instance = auto_ptr<score_matrix> (new score_matrix(daa.score_matrix(),
					daa.gap_open_penalty(),
					daa.gap_extension_penalty(),
					daa.match_reward(),
					daa.mismatch_penalty(),
					_val ()));

	cout << "Build version = " << daa.diamond_build() << endl;
	cout << "DB sequences = " << daa.db_seqs() << endl;
	cout << "DB sequences used = " << daa.db_seqs_used() << endl;
	cout << "DB letters = " << daa.db_letters() << endl;

	View_writer writer;
	const Output_format<_val>& format (get_output_format<_val>());
	format.print_header(writer.f_);

	View_context<_val> context (daa, writer, format);
	// cout<<"View_context..."<<endl;

	launch_thread_pool(context, VATParameters::threads());


}

void view()
{
	DAA_file daa (VATParameters::daa_file);
	if(daa.mode() == blastn)
		view<DNA>(daa);
	else
		view<Protein>(daa);
}

#endif /* VIEW_H_ */
