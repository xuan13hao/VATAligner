

#ifndef LOAD_SEQS_H_
#define LOAD_SEQS_H_

#include <iostream>
#include "SequenceSet.h"
#include "../basic/Translator.h"
#include "../utils/seq_file_format.h"

struct Single_strand { };
struct Double_strand { };

template<typename _ival, typename _val, typename _strand>
size_t push_seq(AlphabetSet<_val> &ss, AlphabetSet<DNA>& source_seqs, const vector<_ival> &seq)
{ ss.push_back(seq); return seq.size(); }

/*
template<>
size_t push_seq<Protein,Nucleotide,Single_strand>(String_set<Nucleotide> &ss, String_set<Nucleotide>& source_seqs, const vector<Protein> &seq)
{ return 0; }

template<>
size_t push_seq<Nucleotide,Protein,Double_strand>(String_set<Protein> &ss, String_set<Nucleotide>& source_seqs, const vector<Nucleotide> &seq)
{
	source_seqs.push_back(seq);
	if(seq.size() < 2) {
		for(unsigned j=0;j<6;++j)
			ss.fill(0, Value_traits<Protein>::MASK_CHAR);
		return 0;
	}
	vector<Protein> proteins[6];
	size_t n = Translator::translate(seq, proteins);

	unsigned bestFrames (Translator::computeGoodFrames(proteins, program_options::get_run_len(seq.size()/3)));
	for(unsigned j = 0; j < 6; ++j) {
		if(bestFrames & (1 << j))
			ss.push_back(proteins[j]);
		else
			ss.fill(proteins[j].size(), Value_traits<Protein>::MASK_CHAR);
	}
	return n;
}

template<>
size_t push_seq<Nucleotide,Nucleotide,Double_strand>(String_set<Nucleotide> &ss, String_set<Nucleotide>& source_seqs, const vector<Nucleotide> &seq)
{
	ss.push_back(seq);
	ss.push_back(Translator::reverse(seq));
	return seq.size()*2;
}*/

template<typename _ival, typename _val, typename _strand>
size_t ReadingSeqs(Input_stream &file,
		const SequenceFileFormat<_ival> &format,
		SequenceSet<_val>** seqs,
		AlphabetSet<char,0>*& ids,
		SequenceSet<DNA>*& source_seqs,
		size_t max_letters)
{
	*seqs = new SequenceSet<_val> ();
	ids = new AlphabetSet<char,0> ();
	source_seqs = new SequenceSet<DNA> ();
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

#endif /* LOAD_SEQS_H_ */