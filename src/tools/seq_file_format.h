
/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef SEQ_FILE_FORMAT_H_
#define SEQ_FILE_FORMAT_H_

using std::pair;

template<typename _val>
class SequenceFileFormat
{
	public:
	virtual bool get_seq(vector<char> &id, vector<_val> &seq, Input_stream &s) const = 0;
	virtual ~SequenceFileFormat()
	{ }

protected:

	template<typename _t>
	static void copy_line(Input_stream &s, vector<_t> &v)
	{
		char a = 0;
		while(s.read(&a, 1) == 1 && a != '\n' && a != '\r')
			v.push_back(AlphabetAttributes<_t>::from_char(a));
		if(a == '\r')
			if(s.read(&a, 1) != 1 || a != '\n')
				throw file_format_exception ();
	}

	static void skip_line(Input_stream &s)
	{
		char a = 0;
		while(s.read(&a, 1) == 1 && a != '\n' && a != '\r');
		if(a == '\r')
			if(s.read(&a, 1) != 1 || a != '\n')
				throw file_format_exception ();
	}

	static void copy_until(Input_stream &s, int delimiter, vector<_val> &v)
	{
		char a = 0;
		size_t col = 0;
		while(s.read(&a, 1) == 1 && a != delimiter) {
			switch(a) {
			case '\n':
				col = 0;
			case '\r':
				break;
			default:
				v.push_back(AlphabetAttributes<_val>::from_char(a));
				++col;
			}
		}
		if(a == delimiter) {
			if(col > 0)
				throw file_format_exception ();
			else
				s.putback(a);
		}
	}

	static void skip_char(Input_stream &s, char c)
	{
		char a;
		if(s.read(&a, 1) != 1 || a != c)
			throw file_format_exception ();
	}

	static bool have_char(Input_stream &s, char c)
	{
		char a;
		if(s.read(&a, 1) == 0)
			return false;
		else if (a != c)
			throw file_format_exception ();
		else
			return true;
	}

};

template<typename _val>
struct FASTA_format : public SequenceFileFormat<_val>
{

	FASTA_format()
	{ }

	virtual bool get_seq(vector<char> &id, vector<_val> &seq, Input_stream &s) const
	{
		if(!SequenceFileFormat<_val>::have_char(s, '>'))
			return false;
		id.clear();
		seq.clear();
		SequenceFileFormat<_val>::copy_line(s, id);
		SequenceFileFormat<_val>::copy_until(s, '>', seq);
		return true;
	}

	virtual ~FASTA_format()
	{ }

};

template<typename _val>
struct FASTQ_format : public SequenceFileFormat<_val>
{

	FASTQ_format()
	{ }

	virtual bool get_seq(vector<char> &id, vector<_val> &seq, Input_stream &s) const
	{
		if(!SequenceFileFormat<_val>::have_char(s, '@'))
			return false;
		id.clear();
		seq.clear();
		SequenceFileFormat<_val>::copy_line(s, id);
		SequenceFileFormat<_val>::copy_line(s, seq);
		SequenceFileFormat<_val>::skip_char(s, '+');
		SequenceFileFormat<_val>::skip_line(s);
		SequenceFileFormat<_val>::skip_line(s);
		return true;
	}

	virtual ~FASTQ_format()
	{ }

};

template<typename _val>
const SequenceFileFormat<_val>* guess_format(const string &file)
{
	static const FASTA_format<_val> fasta;
	static const FASTQ_format<_val> fastq;

	Input_stream f (file, true);
	char c;
	if(f.read(&c, 1) != 1)
		throw file_format_exception ();
	f.close();
	switch(c) {
	case '>': return &fasta;
	case '@': return &fastq;
	default: throw file_format_exception ();
	}
	return 0;
}

#endif /* SEQ_FILE_FORMAT_H_ */
