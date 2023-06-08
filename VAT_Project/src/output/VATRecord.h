

#ifndef DAA_RECORD_H_
#define DAA_RECORD_H_

#include "VATFile.h"
#include "../basic/PackedSequence.h"

using std::string;
using std::vector;

template<typename _val>
struct VATMatchRecord;

template<typename _val>
void translate_query(const vector<DNA>& query, vector<_val> *context)
{
	context[0] = query;
	context[1] = Translator::reverse(query);
}

template<>
void translate_query<Protein>(const vector<DNA>& query, vector<Protein> *context)
{
	Translator::translate(query, context);
}

template<typename _val>
class VATQueryRecord
{
public:
	struct Match_iterator
	{
		Match_iterator(const VATQueryRecord &parent, Binary_buffer::Iterator it):
			r_ (parent),
			it_ (it),
			good_ (true)
		{ operator++(); }
		const VATMatchRecord<_val>& operator*() const
		{ return r_; }
		const VATMatchRecord<_val>* operator->() const
		{ return &r_; }
		bool good() const
		{ return good_; }
		Match_iterator& operator++()
		{ if(it_.good()) it_ >> r_; else good_ = false; return *this; }
	private:
		VATMatchRecord<_val> r_;
		Binary_buffer::Iterator it_;
		bool good_;
	};

	VATQueryRecord(const VATFile& file, const Binary_buffer &buf):
		file_ (file),
		it_ (init(buf))
	{ }

	Match_iterator begin() const
	{ return Match_iterator (*this, it_); }

	string query_name;
	vector<DNA> source_seq;
	vector<_val> context[6];

private:

	Binary_buffer::Iterator init(const Binary_buffer &buf)
	{
		Binary_buffer::Iterator it (buf.begin());
		uint32_t query_len;
		it >> query_len;
		it >> query_name;
		uint8_t flags;
		it >> flags;
		if(file_.mode() == blastp) {
			Packed_sequence seq (it, query_len, false, 5);
			seq.unpack(context[0], 5, query_len);
		} else if(file_.mode() == dna)
		{

			// cout<<"..........."<<endl;
			Packed_sequence seq (it, query_len, false, 5);
			seq.unpack(context[0], 5, query_len);	

						// cout<<"..........."<<endl;
		}
		else {

			const bool have_n = (flags&1) == 1;
			Packed_sequence seq (it, query_len, have_n, have_n ? 3 : 2);
			seq.unpack(source_seq, have_n ? 3 : 2, query_len);
			translate_query<_val>(source_seq, context);
		}
		return it;
	}

	const VATFile& file_;
	const Binary_buffer::Iterator it_;

	friend struct VATMatchRecord<_val>;

	template<typename _val2> friend Binary_buffer::Iterator& operator>>(Binary_buffer::Iterator &it, VATMatchRecord<_val2> &r);

};

template<typename _val>
class VATMatchRecord
{
public:
	VATMatchRecord(const VATQueryRecord<_val> &query_record):
		parent_ (query_record)
	{ }
	bool containsSuffix(const std::string& str, const std::string& suffix) {
		if (str.length() >= suffix.length()) {
			std::string strSuffix = str.substr(str.length() - suffix.length());
			return strSuffix == suffix;
		}
		return false;
	}
	std::string removeSuffix(const std::string& str, const std::string& suffix) {
		if (str.length() >= suffix.length() && str.substr(str.length() - suffix.length()) == suffix) {
			return str.substr(0, str.length() - suffix.length());
		}
		return str;
	}
	const string& query_name() const
	{ return parent_.query_name; }

	const vector<_val>& query() const
	{ return parent_.context[frame]; }

	size_t db_letters() const
	{ return parent_.file_.db_letters(); }

	unsigned query_end() const
	{
		if(parent_.file_.mode() == blastp) {
			return query_begin + translated_query_len - 1;
		} else if(parent_.file_.mode() == blastx) {
			int len = (int)translated_query_len*3*(frame>2 ? -1 : 1);
			return (int)query_begin + (len > 0 ? -1 : 1) + len;
		} else if(parent_.file_.mode() == dna) {
			// int len = (int)translated_query_len*(frame>0 ? -1 : 1);
			// return (int)query_begin + (len > 0 ? -1 : 1) + len;
			return query_begin + translated_query_len - 1;
		} else
			return 0;
	}

	uint32_t total_subject_len, score, query_begin, subject_begin, frame, translated_query_begin, translated_query_len, subject_len, len, identities, mismatches, gap_openings;
	string subject_name;
	Packed_transcript transcript;

private:

	void parse()
	{
		translated_query_len = 0;
		subject_len = 0;
		len = 0;
		identities = 0;
		mismatches = 0;
		gap_openings = 0;
		unsigned d = 0;
		for(Packed_transcript::Const_iterator<char> i = transcript.template begin<char>(); i.good(); ++i) {
			len += i->count;
			switch(i->op) {
			case op_match:
				identities += i->count;
				translated_query_len += i->count;
				subject_len += i->count;
				d = 0;
				break;
			case op_substitution:
				mismatches += i->count;
				translated_query_len += i->count;
				subject_len += i->count;
				d = 0;
				break;
			case op_insertion:
				translated_query_len += i->count;
				++gap_openings;
				d = 0;
				break;
			case op_deletion:
				subject_len += i->count;
				if(d == 0)
					++gap_openings;
				d += i->count;
			}
		}
	}

	const VATQueryRecord<_val> &parent_;

	template<typename _val2> friend Binary_buffer::Iterator& operator>>(Binary_buffer::Iterator &it, VATMatchRecord<_val2> &r);

};

template<typename _val>
Binary_buffer::Iterator& operator>>(Binary_buffer::Iterator &it, VATMatchRecord<_val> &r)
{
	uint32_t subject_id;
	it >> subject_id;
	uint8_t flag;
	it >> flag;
	it.read_packed(flag & 3, r.score);
	it.read_packed((flag>>2)&3, r.query_begin);
	it.read_packed((flag>>4)&3, r.subject_begin);
	r.transcript.read(it);
	r.subject_name = r.parent_.file_.ref_name(subject_id);
	r.total_subject_len = r.parent_.file_.ref_len(subject_id);
	// cout<<"r.subject_name = "<<r.subject_name<<",r.total_subject_len =  "<<r.total_subject_len<<endl;
	if(r.parent_.file_.mode() == blastx) {
		r.frame = (flag&(1<<6)) == 0 ? r.query_begin % 3 : 3+(r.parent_.source_seq.size() - 1 - r.query_begin)%3;
		r.translated_query_begin = query_translated_begin<_val>(r.query_begin, r.frame, r.parent_.source_seq.size(), true);
	} else if (r.parent_.file_.mode() == blastp) {
		r.frame = 0;
		r.translated_query_begin = r.query_begin;
	} else {
		r.frame = 0;
		r.translated_query_begin = r.query_begin;
		// if(r.containsSuffix(r.subject_name ,"_minus"))
		// {
		// 	r.subject_begin = r.total_subject_len - r.len;
		// }
		// r.frame = (flag&(1<<6)) == 0 ? 0 : 1;
		// r.translated_query_begin = query_translated_begin<_val>(r.query_begin, r.frame, r.parent_.source_seq.size(), false);
	}
	r.parse();
	return it;
}




#endif /* DAA_RECORD_H_ */
