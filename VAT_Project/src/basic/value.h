

#ifndef VALUE_H_
#define VALUE_H_

#include "AlphabetType.h"
#include "VATConsts.h"

typedef enum { amino_acid=0, nucleotide=1 } Sequence_type;

template<typename _val>
Sequence_type sequence_type(const _val&)
{ return amino_acid; }

template<>
Sequence_type sequence_type<DNA>(const DNA&)
{ return nucleotide; }

Sequence_type input_sequence_type()
{ return VATParameters::algn_type == VATParameters::protein ? amino_acid : nucleotide; }

Sequence_type sequence_type()
{ return VATParameters::algn_type == VATParameters::dna ? nucleotide : amino_acid; }

size_t query_contexts()
{
	// switch(VATParameters::command) {
	// case VATParameters::dna: return 1;
	// case VATParameters::blastx: return 6;
	// default: return 1;
	// }
	return 1;
}

bool query_translated()
{ return VATParameters::algn_type == VATParameters::blastx ? true : false; }

int query_len_factor()
{ return VATParameters::algn_type == VATParameters::blastx ? 3 : 1; }

template<typename _val>
class AlphabetMap
{
	public:
	AlphabetMap(unsigned size, const char *chars, char mask, const char *mask_chars)
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

template<> const Protein AlphabetMap<Protein>::invalid = 0xff;
template<> const DNA AlphabetMap<DNA>::invalid = 0xff;

template<typename _val>
struct AlphabetAttributes
{ };

template<>
struct AlphabetAttributes<Protein>
{
	enum { ALPHABET_SIZE = 25 };
	static const Protein				MASK_CHAR;
	static const char*					ALPHABET;
	static const AlphabetMap<Protein>	from_char;
};

const Protein					AlphabetAttributes<Protein>::MASK_CHAR = 23;
const char* AlphabetAttributes<Protein>::ALPHABET = "ARNDCQEGHILKMFPSTWYVBJZX*";
const AlphabetMap<Protein> AlphabetAttributes<Protein>::from_char (AlphabetAttributes<Protein>::ALPHABET_SIZE, AlphabetAttributes<Protein>::ALPHABET, AlphabetAttributes<Protein>::MASK_CHAR, "UO-");

template<>
struct AlphabetAttributes<const Protein> : public AlphabetAttributes<Protein>
{ };

template<>
struct AlphabetAttributes<DNA>
{
	enum { ALPHABET_SIZE = 5 };
	static const DNA				MASK_CHAR;
	static const char*					ALPHABET;
	static const AlphabetMap<DNA>	from_char;
};

const DNA AlphabetAttributes<DNA>::MASK_CHAR = 4;
const char* AlphabetAttributes<DNA>::ALPHABET = "ACGTN";
const AlphabetMap<DNA> AlphabetAttributes<DNA>::from_char (AlphabetAttributes<DNA>::ALPHABET_SIZE, AlphabetAttributes<DNA>::ALPHABET, AlphabetAttributes<DNA>::MASK_CHAR, "MRWSYKVHDBX");

template<>
struct AlphabetAttributes<const DNA> : public AlphabetAttributes<DNA>
{ };

char to_char(Protein a)
{ return AlphabetAttributes<Protein>::ALPHABET[a]; }

template<>
struct AlphabetAttributes<char>
{
	static char from_char(char c)
	{ return c; }
};

#endif /* VALUE_H_ */
