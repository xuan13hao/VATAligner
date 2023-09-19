#ifndef __ALPHABETFEATURE_H__
#define __ALPHABETFEATURE_H__
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include "SequenceType.h"
#include "BioAlphabet.h"
using std::cout;
using std::endl;


template<typename T>
class AlphabetFeature
{
   
};

template<>
class AlphabetFeature<DNA>
{
    public:
    /* data */
    enum
    {
        ALPHABET_SIZE = 5
    };
    static const DNA mask_char;//MASK_CHAR
    static const char *alpha;//ALPHABET
    static const BioAlphabet<DNA> bioa;//type of sequence 
};

template<>
class AlphabetFeature<RNA>
{
    public:
    /* data */
    enum
    {
        ALPHABET_SIZE = 5
    };
    static const RNA mask_char;//MASK_CHAR
    static const char *alpha;//ALPHABET
    static const BioAlphabet<RNA> bioa;//type of sequence 
};

template<>
class AlphabetFeature<Protein>
{
    public:
    
	enum 
    { 
        ALPHABET_SIZE = 25 
    };
	static const Protein				mask_char;
	static const char*					alpha;
	static const BioAlphabet<Protein>	bioa;
};

template<>
class AlphabetFeature<const DNA>: public AlphabetFeature<DNA>
{

};

template<>
class AlphabetFeature<const RNA>: public AlphabetFeature<RNA>
{

};

template<>
class AlphabetFeature<const Protein> : public AlphabetFeature<Protein>
{ 

};

template<>
class AlphabetFeature<char>
{
    public:
    static char bioa(char str)
    {
        return str;
    }
};


#endif // __ALPHABETFEATURE_H__