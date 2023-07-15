
#ifndef OUTPUT_FORMAT_H_
#define OUTPUT_FORMAT_H_

#include "../basic/Hits.h"
#include "../align/match_func.h"
#include "../output/VATFile.h"
#include "../output/VATRecord.h"

template<typename _val>
class Output_format
{
	public:
	virtual void print_match(const VATMatchRecord<_val> &r, Text_buffer &out) const = 0;
	virtual void print_header(OutputStreamer &f) const
	{ }
	virtual ~Output_format()
	{ }

};

template<typename _val>
class Blast_tab_format : public Output_format<_val>
{
	public:
	Blast_tab_format()
	{ }
	virtual void print_match(const VATMatchRecord<DNA> &r, Text_buffer &out) const
	{



		if ((containsSuffix(r.query_name(),"_minus")))
		{
			string qry_name = removeSuffix(r.query_name(),"_minus");
			uint32_t sbj_start= r.subject_begin+r.subject_len;
			uint32_t sbj_end= r.subject_begin+1;
			// cout<<"sub start = "<<r.subject_begin<<",r.subject_len =  "<<r.subject_len<<", r.total_subject_len = "<<r.total_subject_len<<endl;
			// cout<<"sub start = "<<sbj_start<<",sub end =  "<<sbj_end<<", r.total_subject_len = "<<r.total_subject_len<<endl;
			out << qry_name << '\t'
					<< r.subject_name << '\t'
					<< (double)r.identities*100/r.len << '\t'
					<< r.len << '\t'
					<< r.mismatches << '\t'
					<< r.gap_openings << '\t'
					// << "Plus/Minus" << '\t'
					<< r.query_begin+1 << '\t'
					<< r.query_end()+1 << '\t'
					<< sbj_start<< '\t'
					<< sbj_end<< '\t';
			out.print_e(ScoreMatrix::get().evalue(r.score, r.db_letters(), r.query().size()));
			out << '\t' << ScoreMatrix::get().bitscore(r.score) << '\n';
		}
		else{		
			out << r.query_name() << '\t'
					<< r.subject_name << '\t'
					<< (double)r.identities*100/r.len << '\t'
					<< r.len << '\t'
					<< r.mismatches << '\t'
					<< r.gap_openings << '\t'
					// << "Plus/Plus" << '\t'
					<< r.query_begin+1 << '\t'
					<< r.query_end()+1 << '\t'
					<< r.subject_begin+1<< '\t'
					<< r.subject_begin+r.subject_len << '\t';
			out.print_e(ScoreMatrix::get().evalue(r.score, r.db_letters(), r.query().size()));
			out << '\t' << ScoreMatrix::get().bitscore(r.score) << '\n';
		}
	}

		virtual void print_match(const VATMatchRecord<Protein> &r, Text_buffer &out) const
	{
			out << r.query_name() << '\t'
					<< r.subject_name << '\t'
					<< (double)r.identities*100/r.len << '\t'
					<< r.len << '\t'
					<< r.mismatches << '\t'
					<< r.gap_openings << '\t'
					<< r.query_begin+1 << '\t'
					<< r.query_end()+1 << '\t'
					<< r.subject_begin+1<< '\t'
					<< r.subject_begin+r.subject_len << '\t';
			out.print_e(ScoreMatrix::get().evalue(r.score, r.db_letters(), r.query().size()));
			out << '\t' << ScoreMatrix::get().bitscore(r.score) << '\n';
	}

	virtual ~Blast_tab_format()
	{ }
	static bool containsSuffix(const std::string& str, const std::string& suffix) {
		if (str.length() >= suffix.length()) {
			std::string strSuffix = str.substr(str.length() - suffix.length());
			return strSuffix == suffix;
		}
		return false;
	}
	static std::string removeSuffix(const std::string& str, const std::string& suffix) {
		if (str.length() >= suffix.length() && str.substr(str.length() - suffix.length()) == suffix) {
			return str.substr(0, str.length() - suffix.length());
		}
		return str;
	}
	static size_t print_salltitles(Text_buffer &buf, const char *id)
	{
		size_t n = 0;
		const vector<string> t (tokenize(id, "\1"));
		vector<string>::const_iterator i=t.begin();
		for(;i<t.end()-1;++i) {
			buf << *i << "<>";
			n += i->length() + 2;
		}
		buf << *i;
		n += i->length();
		return n;
	}

};

template<typename _val>
class Sam_format : public Output_format<_val>
{
	public:
	Sam_format()
	{ }

	virtual void print_match(const VATMatchRecord<_val> &r, Text_buffer &out) const
	{
		out << r.query_name() << '\t'
				<< '0' << '\t'
				<< r.subject_name << '\t'
				<< r.subject_begin+1 << '\t'
				<< "255" << '\t';

		print_cigar(r, out);

		out << '\t'
				<< '*' << '\t'
				<< '0' << '\t'
				<< '0' << '\t'
				<< sequence<const _val> (&r.query()[r.translated_query_begin], r.translated_query_len) << '\t'
				<< '*' << '\t'
				<< "AS:i:" << (uint32_t)ScoreMatrix::get().bitscore(r.score) << '\t'
				<< "NM:i:" << r.len - r.identities << '\t'
				<< "ZL:i:" << r.total_subject_len << '\t'
				<< "ZR:i:" << r.score << '\t'
				<< "ZE:f:";
		out.print_e(ScoreMatrix::get().evalue(r.score, r.db_letters(), r.query().size()));
		out << '\t'
				<< "ZI:i:" << r.identities*100/r.len << '\t'
				<< "ZF:i:" << blast_frame(r.frame) << '\t'
				<< "ZS:i:" << r.query_begin+1 << '\t'
				<< "MD:Z:";

		print_md(r, out);
		out << '\n';
	}

	void print_md(const VATMatchRecord<_val> &r, Text_buffer &buf) const
	{
		unsigned matches = 0, del = 0;
		for(Packed_transcript::Const_iterator<_val> i = r.transcript.template begin<_val>(); i.good(); ++i) {
			switch(i->op) {
			case op_match:
				del = 0;
				matches += i->count;
				break;
			case op_insertion:
				break;
			case op_substitution:
				if(matches > 0) {
					buf << matches;
					matches = 0;
				} else if(del > 0) {
					buf << '0';
					del = 0;
				}
				buf << AlphabetAttributes<_val>::ALPHABET[i->letter];
				break;
			case op_deletion:
				if(matches > 0) {
					buf << matches;
					matches = 0;
				}
				if(del == 0)
					buf << '^';
				buf << AlphabetAttributes<_val>::ALPHABET[i->letter];
				++del;
			}
		}
		if(matches > 0)
			buf << matches;
	}

	void print_cigar(const VATMatchRecord<_val> &r, Text_buffer &buf) const
	{
		static const unsigned map[] = { 0, 1, 2, 0 };
		static const char letter[] = { 'M', 'I', 'D' };
		unsigned n = 0, op = 0;
		for(Packed_transcript::Const_iterator<_val> i = r.transcript.template begin<_val>(); i.good(); ++i) {
			if(map[i->op] == op)
				n += i->count;
			else {
				if(n > 0)
					buf << n << letter[op];
				n = i->count;
				op = map[i->op];
			}
		}
		if(n > 0)
			buf << n << letter[op];
	}

	virtual void print_header(OutputStreamer &f) const
	{
		static const char* line = "@HD\tVN:1.5\tSO:query\n\
@PG\tPN:VAT\n\
@mm\tVAT\n\
@CO\tVAT alignments\n\
@CO\tReporting AS: bitScore, ZR: rawScore, ZE: expected, ZI: percent identity, ZL: reference length, ZF: frame, ZS: query start DNA coordinate\n";
		f.write(line, strlen(line));
	}

	virtual ~Sam_format()
	{ }

};

template<typename _val>
const Output_format<_val>& get_output_format()
{
	static const Sam_format<_val> sam;
	static const Blast_tab_format<_val> tab;
	if(VATParameters::output_format == "tab")
		return tab;
	else if(VATParameters::output_format == "sam")
		return sam;
	else
		throw std::runtime_error("Invalid output format.");
}

#endif /* OUTPUT_FORMAT_H_ */
