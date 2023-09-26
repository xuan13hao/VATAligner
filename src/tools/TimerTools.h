
#ifndef LOG_STREAM_H_
#define LOG_STREAM_H_

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/timer/timer.hpp>

using std::string;
using std::endl;
using std::cout;
boost::iostreams::filtering_ostream verbose_stream;
boost::iostreams::filtering_ostream log_stream;

class TimerTools : public boost::timer::cpu_timer
{
	public:
	TimerTools(const char *msg, bool print):
		print_ (print),
		msg_ (msg)
	{ start(msg); }
	~TimerTools()
	{ finish(); }
	void go(const char *msg)
	{
		finish();
		boost::timer::cpu_timer::start();
		start(msg);
		msg_ = msg;
	}
	void finish()
	{
		if(!msg_)
			return;
		if(print_ && !VATParameters::debug_log)
			log_stream << boost::timer::format(elapsed(), 1, "[%ws]") << endl;
		else {
			verbose_stream << ' ' << msg_ << boost::timer::format(elapsed(), 1, " [%ws]") << endl;
		}
		msg_ = 0;
	}
private:
	void start(const char *msg)
	{
		if(print_ && !VATParameters::debug_log) {
			log_stream << msg << "... " << std::flush;
			fflush(stdout);
		} else
			verbose_stream << msg << "..." << endl;
	}
	bool print_;
	const char *msg_;
};

#endif
