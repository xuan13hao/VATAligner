#ifndef __BASEMATRIX_H__
#define __BASEMATRIX_H__

#include <math.h>


//template <typename T>
class NucleotideMatrix 
{
    private:
        const int  a =  1;   // Match
        const int  b =  -4;   // Mismatch
        const int gap_open = -5;
        const int gap_ext = -3;
        const int blast_lambbda = 1.28;
        const double blast_k = 0.46;


    public:
    /**
     * 
     * double bit_score = (NT_BLAST_LAMBDA * result[i][j].score - log(NT_BLAST_K)) / log(2);
       double evalue = log10(strlen(seq[begin + i])) + log10(db_size) - bit_score * log10(2);
       int e_exponent = floor(evalue);
    */


        double evalue(double bit_score, int db_letters, unsigned query_len) const
        {
            double evalue = log10(query_len) + log10(db_letters) - bit_score * log10(2);
            return floor(evalue);
        }
        double bitscore(double score) const
        {
            double bit_score = (blast_lambbda * score - log(blast_k)) / log(2);
            return bit_score;
        }

};


#endif // __BASEMATRIX_H__