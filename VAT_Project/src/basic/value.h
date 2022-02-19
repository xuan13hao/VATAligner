

#ifndef VALUE_H_
#define VALUE_H_

#include "value_type.h"
#include "VATConsts.h"

typedef enum { amino_acid=0, nucleotide=1 } Sequence_type;

template<typename _val>
Sequence_type sequence_type(const _val&)
{ return amino_acid; }

template<>
Sequence_type sequence_type<DNA>(const DNA&)
{ return nucleotide; }

Sequence_type input_sequence_type()
{ return VATParameters::command == VATParameters::blastp ? amino_acid : nucleotide; }

Sequence_type sequence_type()
{ return VATParameters::command == VATParameters::blastn ? nucleotide : amino_acid; }

size_t query_contexts()
{
	switch(VATParameters::command) {
	case VATParameters::blastn: return 2;
	case VATParameters::blastx: return 6;
	default: return 1;
	}
}

bool query_translated()
{ return VATParameters::command == VATParameters::blastx ? true : false; }

int query_len_factor()
{ return VATParameters::command == VATParameters::blastx ? 3 : 1; }

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

template<> const Protein Char_representation<Protein>::invalid = 0xff;
template<> const DNA Char_representation<DNA>::invalid = 0xff;

template<typename _val>
struct Value_traits
{ };

template<>
struct Value_traits<Protein>
{
	enum { ALPHABET_SIZE = 25 };
	static const Protein				MASK_CHAR;
	static const char*					ALPHABET;
	static const Char_representation<Protein>	from_char;
};

const Protein					Value_traits<Protein>::MASK_CHAR = 23;
const char* Value_traits<Protein>::ALPHABET = "ARNDCQEGHILKMFPSTWYVBJZX*";
const Char_representation<Protein> Value_traits<Protein>::from_char (Value_traits<Protein>::ALPHABET_SIZE, Value_traits<Protein>::ALPHABET, Value_traits<Protein>::MASK_CHAR, "UO-");

template<>
struct Value_traits<const Protein> : public Value_traits<Protein>
{ };

template<>
struct Value_traits<DNA>
{
	enum { ALPHABET_SIZE = 5 };
	static const DNA				MASK_CHAR;
	static const char*					ALPHABET;
	static const Char_representation<DNA>	from_char;
};

const DNA Value_traits<DNA>::MASK_CHAR = 4;
const char* Value_traits<DNA>::ALPHABET = "ACGTN";
const Char_representation<DNA> Value_traits<DNA>::from_char (Value_traits<DNA>::ALPHABET_SIZE, Value_traits<DNA>::ALPHABET, Value_traits<DNA>::MASK_CHAR, "MRWSYKVHDBX");

template<>
struct Value_traits<const DNA> : public Value_traits<DNA>
{ };

char to_char(Protein a)
{ return Value_traits<Protein>::ALPHABET[a]; }

template<>
struct Value_traits<char>
{
	static char from_char(char c)
	{ return c; }
};

#endif /* VALUE_H_ */
