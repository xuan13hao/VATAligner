
#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include <string>
#include <boost/lexical_cast.hpp>
#include "../tools/tinythread.h"

using std::exception;
using std::string;

#define THROW_EXCEPTION(exception, p1) throw exception(p1, __PRETTY_FUNCTION__, __LINE__)

struct VATException : public exception
{
	string msg;
	VATException(const char *function, unsigned line):
		msg(string("function ") + function + " line " + boost::lexical_cast<string>(line))
	{ }
	VATException(const string &msg):
		msg (msg)
	{ }
	~VATException() throw() { }
	virtual const char* what() const throw()
	{ return msg.c_str(); }
};

struct file_io_exception : public VATException
{
	file_io_exception(const string &file_name, const char* function, unsigned line):
		VATException(function, line)
	{ msg += ". Error reading file " + file_name; }
};

struct file_open_exception: public VATException
{
	file_open_exception(const string &file_name, const char* function, unsigned line):
		VATException(function, line)
	{ msg += ". Error opening file " + file_name; }
};

struct file_io_write_exception: public VATException
{
	file_io_write_exception(const string &file_name, const char* function, unsigned line):
		VATException(function, line)
	{ msg += ". Error writing file " + file_name; }
};

struct memory_alloc_exception: public exception
{
	virtual const char* what() const throw()
	{ return "Failed to allocate memory"; }
};

struct hash_table_overflow_exception : public exception
{
	virtual const char* what() const throw()
	{ return "Hash table overflow"; }
};

struct invalid_database_version_exception : public exception
{
	virtual const char* what() const throw()
	{ return "Incompatible database version"; }
};

struct invalid_parameter_exception : public exception
{
	virtual const char* what() const throw()
	{ return "Invalid parameter"; }
};

struct invalid_sequence_char_exception : public exception
{
	const string msg;
	invalid_sequence_char_exception(char ch):
		msg (string("Invalid character (")+ch+"/" + boost::lexical_cast<string>(int(ch)) + ") in sequence")
	{ }
	~invalid_sequence_char_exception() throw()
	{ }
	virtual const char* what() const throw()
	{ return msg.c_str(); }
};

struct file_format_exception : public exception
{
	virtual const char* what() const throw()
	{ return "Invalid input file format"; }
};

struct Exception_state
{
	Exception_state():
		active_ (false)
	{ }
	void set(const std::exception &e)
	{
		if(!active_) {
			mtx_.lock();
			what_ = string(typeid(e).name()) + ": " + e.what();
			active_ = true;
			mtx_.unlock();
		}
	}
	void sync() const
	{
		if(active_)
			throw std::runtime_error(what_);
	}
	bool operator()() const
	{ return active_; }
private:
	tthread::mutex mtx_;
	bool active_;
	string what_;
} exception_state;

#endif /* EXCEPTIONS_H_ */
