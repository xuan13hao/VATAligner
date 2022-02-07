

#ifndef LOAD_SEQS_H_
#define LOAD_SEQS_H_

#include <iostream>
#include "sequence_set.h"
#include "../basic/translate.h"
#include "../util/seq_file_format.h"

struct Single_strand { };
struct Double_strand { };
//<DNA,_val,Single_strand>

//ival = dna;
template <typename _ival, typename _val, typename _strand>
size_t push_seq(String_set<_val> &ss, String_set<DNA> &source_seqs, const vector<_ival> &seq)
{
	//cout << "1" << endl;
	ss.push_back(seq);
	return seq.size();
}

template<typename _ival, typename _val, typename _strand>
size_t push_seq(String_set<_val> &ss, String_set<RNA>& source_seqs, const vector<_ival> &seq)
{

	ss.push_back(seq);
	return seq.size();
}


size_t pushDNASeq(String_set<DNA> &ss, String_set<DNA>& source_seqs, const vector<DNA> &seq)
{
	ss.push_back(seq);
	return seq.size();
}


size_t pushRNASeq(String_set<RNA> &ss, String_set<RNA>& source_seqs, const vector<RNA> &seq)
{
	ss.push_back(seq);
	return seq.size();
}


size_t pushProteinSeq(String_set<Protein> &ss, String_set<DNA>& source_seqs, const vector<Protein> &seq)
{
	ss.push_back(seq);
	return seq.size();
}

// template<>
// size_t push_seq<Protein,DNA,Single_strand>(String_set<DNA> &ss, String_set<DNA>& source_seqs, const vector<Protein> &seq)
// {
// 	cout << "push seq 2" << endl;
// 	return 0;
// }

template<>
size_t push_seq<DNA,Protein,Double_strand>(String_set<Protein> &ss, String_set<DNA>& source_seqs, const vector<DNA> &seq)
{
	source_seqs.push_back(seq);
	if(seq.size() < 2) {
		for(unsigned j=0;j<6;++j)
			ss.fill(0, AlphabetFeature<Protein>::mask_char);
		return 0;
	}
	vector<Protein> proteins[6];
	size_t n = Translator::translate(seq, proteins);

	unsigned bestFrames (Translator::computeGoodFrames(proteins, program_options::get_run_len(seq.size()/3)));
	for(unsigned j = 0; j < 6; ++j) {
		if(bestFrames & (1 << j))
			ss.push_back(proteins[j]);
		else
			ss.fill(proteins[j].size(), AlphabetFeature<Protein>::mask_char);
	}
	return n;
}

template<>
size_t push_seq<DNA,DNA,Double_strand>(String_set<DNA> &ss, String_set<DNA>& source_seqs, const vector<DNA> &seq)
{
	ss.push_back(seq);
	ss.push_back(Translator::reverse(seq));
	return seq.size()*2;
}

template<typename _ival, typename _val, typename _strand>
size_t load_seqs(Input_stream &file,
		const Sequence_file_format<_ival> &format,
		Sequence_set<_val>** seqs,
		String_set<char,0>*& ids,
		Sequence_set<DNA>*& source_seqs,
		size_t max_letters)
{
	*seqs = new Sequence_set<_val> ();
	ids = new String_set<char,0> ();
	source_seqs = new Sequence_set<DNA> ();
	size_t letters = 0, n = 0;
	vector<_ival> seq;
	vector<char> id;
	try {
		while(letters < max_letters && format.get_seq(id, seq, file)) {
			ids->push_back(id);
			letters += push_seq<_ival,_val,_strand>(**seqs, *source_seqs, seq);
			++n;
		}
	} catch(invalid_sequence_char_exception &e) {
		std::cerr << n << endl;
		throw e;
	}
	ids->finish_reserve();
	(*seqs)->finish_reserve();
	source_seqs->finish_reserve();
	if(n == 0) {
		delete *seqs;
		delete ids;
		delete source_seqs;
	}
	return n;
}



template<typename _ival, typename _val, typename _strand>
size_t loadDNASeqs(Input_stream &file,
		const Sequence_file_format<_ival> &format,
		Sequence_set<_val>** seqs,
		String_set<char,0>*& ids,
		Sequence_set<DNA>*& source_seqs,
		size_t max_letters)
{
	*seqs = new Sequence_set<_val> ();
	ids = new String_set<char,0> ();
	source_seqs = new Sequence_set<DNA> ();
	size_t letters = 0, n = 0;
	vector<_ival> seq;
	vector<char> id;
	try {
		while(letters < max_letters && format.get_seq(id, seq, file)) {
			ids->push_back(id);
			letters += pushDNASeq(**seqs, *source_seqs, seq);
			++n;
		}
	} catch(invalid_sequence_char_exception &e) {
		std::cerr << n << endl;
		throw e;
	}
	ids->finish_reserve();
	(*seqs)->finish_reserve();
	source_seqs->finish_reserve();
	if(n == 0) {
		delete *seqs;
		delete ids;
		delete source_seqs;
	}
	return n;
}

template<typename _ival, typename _val, typename _strand>
size_t loadProteinSeqs(Input_stream &file,
		const Sequence_file_format<_ival> &format,
		Sequence_set<_val>** seqs,
		String_set<char,0>*& ids,
		Sequence_set<DNA>*& source_seqs,
		size_t max_letters)
{
	*seqs = new Sequence_set<_val> ();
	ids = new String_set<char,0> ();
	source_seqs = new Sequence_set<DNA> ();
	size_t letters = 0, n = 0;
	vector<_ival> seq;
	vector<char> id;
	try {
		while(letters < max_letters && format.get_seq(id, seq, file)) {
			ids->push_back(id);
			letters += pushProteinSeq(**seqs, *source_seqs, seq);
			++n;
		}
	} catch(invalid_sequence_char_exception &e) {
		std::cerr << n << endl;
		throw e;
	}
	ids->finish_reserve();
	(*seqs)->finish_reserve();
	source_seqs->finish_reserve();
	if(n == 0) {
		delete *seqs;
		delete ids;
		delete source_seqs;
	}
	return n;
}

template<typename _ival, typename _val, typename _strand>
size_t loadRNAseqs(Input_stream &file,
		const Sequence_file_format<_ival> &format,
		Sequence_set<_val>** seqs,
		String_set<char,0>*& ids,
		Sequence_set<RNA>*& source_seqs,
		size_t max_letters)
{
	*seqs = new Sequence_set<_val> ();
	ids = new String_set<char,0> ();
	source_seqs = new Sequence_set<RNA> ();
	size_t letters = 0, n = 0;
	vector<_ival> seq;
	vector<char> id;
	try {
		while(letters < max_letters && format.get_seq(id, seq, file)) {
			ids->push_back(id);
			letters += pushRNASeq(**seqs, *source_seqs, seq);
			++n;
		}
	} catch(invalid_sequence_char_exception &e) {
		std::cerr << n << endl;
		throw e;
	}
	ids->finish_reserve();
	(*seqs)->finish_reserve();
	source_seqs->finish_reserve();
	if(n == 0) {
		delete *seqs;
		delete ids;
		delete source_seqs;
	}
	return n;
}

#endif /* LOAD_SEQS_H_ */
