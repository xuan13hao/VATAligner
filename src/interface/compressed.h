#ifndef __COMPRESSED_H__
#define __COMPRESSED_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "../basic/exceptions.h"
#include "temp_file.h"

namespace io = boost::iostreams;

#include <iostream>
#include <fstream>
#include <string>
struct OutputStreamer : public io::filtering_ostream
{

	enum Flags { file_sink, stdout_sink };

	OutputStreamer()
	{ }

	OutputStreamer(const string &file_name, bool gzipped = false, Flags flags = file_sink):
		file_name_ (file_name)
	{
		if(gzipped)
			this->push(io::gzip_compressor ());
		if(flags == file_sink) {
			io::file_sink f (file_name, std::ios_base::out | std::ios_base::binary);
			if(!f.is_open())
				THROW_EXCEPTION(file_open_exception, file_name_);
			this->push(f);
		} else
			this->push(std::cout);
	}

	OutputStreamer(const TempFile &tmp_file)
	{
		this->push(io::file_descriptor_sink(tmp_file.file_descriptor(), io::never_close_handle));
	}

	template<typename _t>
	void write(const _t *ptr, size_t count)
	{
		ssize_t n;
		if((n = io::write(*this, reinterpret_cast<const char*>(ptr), sizeof(_t) * count)) != (ssize_t)(sizeof(_t) * count))
			throw File_write_exception (file_name_.c_str(), sizeof(_t) * count, n);
	}

	void write_c_str(const string &s)
	{ write(s.c_str(), s.length()+1); }

	template<class _t>
	void write(const vector<_t> &v, bool write_size = true)
	{
		size_t size = v.size();
		if(write_size)
			write(&size, 1);
		write(v.data(), size);
	}

	size_t tell()
	{
		if(this->tellp() < 0)
			throw std::runtime_error("Unable to execute tellp() on output stream.");
		return this->tellp();
	}

	void seek(size_t pos)
	{
		this->seekp(pos);
		if(!this->good())
			throw std::runtime_error("Unable to execute seekp() on output stream.");
	}

	void close()
	{
		this->set_auto_close(true);
		this->pop();
	}

private:

	const string file_name_;

};

int main() {
    // Example usage
    std::string data = "Example data";

    OutputStreamer output(true);  // Use compression

    output.write(data);

    return 0;
}

#endif // __COMPRESSED_H__