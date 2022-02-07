
#ifndef VALUE_H_
#define VALUE_H_

// #include "value_type.h"
#include "AlphabetFeature.h"
#include "BioAlphabet.h"
#include "VATParameter.h"


typedef enum { amino_acid=0, nucleotide=1 } Sequence_type;

template<typename _val>
Sequence_type sequence_type(const _val&)
{ return amino_acid; }

template<>
Sequence_type sequence_type<DNA>(const DNA&)
{ return nucleotide; }

Sequence_type input_sequence_type()
{ return program_options::command == program_options::blastp ? amino_acid : nucleotide; }

Sequence_type sequence_type()
{ return program_options::command == program_options::blastn ? nucleotide : amino_acid; }

size_t query_contexts()
{
	switch(program_options::command) {
	case program_options::blastn: return 2;
	case program_options::blastx: return 6;
	default: return 1;
	}
}

bool query_translated()
{ return program_options::command == program_options::blastx ? true : false; }

int query_len_factor()
{ return program_options::command == program_options::blastx ? 3 : 1; }


template <> const DNA BioAlphabet<DNA>::invalid = 0xff;
template <> const Protein BioAlphabet<Protein>::invalid = 0xff;
template <> const RNA BioAlphabet<RNA>::invalid = 0xff;


const DNA AlphabetFeature<DNA>::mask_char = 4;
const char *AlphabetFeature<DNA>::alpha = "ACGTN";
const Protein AlphabetFeature<Protein>::mask_char = 23;
const char* AlphabetFeature<Protein>::alpha = "ARNDCQEGHILKMFPSTWYVBJZX*";
const BioAlphabet<DNA> AlphabetFeature<DNA>::bioa('D');
const BioAlphabet<Protein> AlphabetFeature<Protein>::bioa('P');


/*
template<typename _val>
struct Char_representation
{
	Char_representation(unsigned size, const char *chars, char mask, const char *mask_chars)
	{
		memset(data_, invalid, sizeof(data_));
		for(unsigned i=0;i<size;++i) {
			assert(chars[i] != (char)invalid);
			data_[(long)chars[i]] = i;
			data_[(long)tolower(chars[i])] = i;
		}
		while(*mask_chars != 0) {
			const char ch = *mask_chars;
			data_[(long)ch] = mask;
			data_[(long)tolower(ch)] = mask;
			++mask_chars;
		}
	}
	_val operator()(char c) const
	{
		if(data_[(long)c] == invalid)
			throw invalid_sequence_char_exception (c);
		return data_[(long)c];
	}
private:
	static const _val invalid;
	_val data_[256];
};

template<> const Amino_acid Char_representation<Amino_acid>::invalid = 0xff;
template<> const Nucleotide Char_representation<Nucleotide>::invalid = 0xff;

template<typename _val>
struct Value_traits
{ };

template<>
struct Value_traits<Amino_acid>
{
	enum { ALPHABET_SIZE = 25 };
	static const Amino_acid				MASK_CHAR;
	static const char*					ALPHABET;
	static const Char_representation<Amino_acid>	from_char;
};

const Amino_acid					Value_traits<Amino_acid>::MASK_CHAR = 23;
const char* Value_traits<Amino_acid>::ALPHABET = "ARNDCQEGHILKMFPSTWYVBJZX*";
const Char_representation<Amino_acid> Value_traits<Amino_acid>::from_char (Value_traits<Amino_acid>::ALPHABET_SIZE, Value_traits<Amino_acid>::ALPHABET, Value_traits<Amino_acid>::MASK_CHAR, "UO-");

template<>
struct Value_traits<const Amino_acid> : public Value_traits<Amino_acid>
{ };

template<>
struct Value_traits<Nucleotide>
{
	enum { ALPHABET_SIZE = 5 };
	static const Nucleotide				MASK_CHAR;
	static const char*					ALPHABET;
	static const Char_representation<Nucleotide>	from_char;
};

const Nucleotide Value_traits<Nucleotide>::MASK_CHAR = 4;
const char* Value_traits<Nucleotide>::ALPHABET = "ACGTN";
const Char_representation<Nucleotide> Value_traits<Nucleotide>::from_char (Value_traits<Nucleotide>::ALPHABET_SIZE, Value_traits<Nucleotide>::ALPHABET, Value_traits<Nucleotide>::MASK_CHAR, "MRWSYKVHDBX");

template<>
struct Value_traits<const Nucleotide> : public Value_traits<Nucleotide>
{ };

char to_char(Amino_acid a)
{ return Value_traits<Amino_acid>::ALPHABET[a]; }

template<>
struct Value_traits<char>
{
	static char from_char(char c)
	{ return c; }
};
*/
#endif /* VALUE_H_ */
