#ifndef __NUCLSCOREMATRIX_H__
#define __NUCLSCOREMATRIX_H__


#include "../algo/blast/core/blast_stat.h"
#include "../algo/blast/core/blast_encoding.h"
#include "ProteinProfile.h"
// #include "Matrixs.h"
#include "AlphabetType.h"
#include <memory>

using std::string;
using std::cout;
using std::endl;
using std::auto_ptr;



class ScoreParamsException : public std::exception
{
	public:
	virtual const char* what() const throw()
	{ 
		return "Error Scoring Parameters"; 
	}
};

class Blastscoreblk
{
	public:
	template<typename _val>
	Blastscoreblk(const string &matrix, int gap_open, int gap_extend, int reward, int penalty, const _val&):
		data_ (BlastScoreBlkNew(blast_seq_code<_val>(), 1))
	{
		if(data_ == 0)
			throw ScoreParamsException ();
		if((data_->kbp_gap_std[0] = Blast_KarlinBlkNew()) == 0
				|| (data_->kbp_std[0] = Blast_KarlinBlkNew()) == 0)
			throw ScoreParamsException ();
		if(blast_load_karlin_blk<_val>(data_->kbp_gap_std[0],
				data_->kbp_std[0],
				gap_open,
				gap_extend,
				reward,
				penalty,
				matrix.c_str()) != 0)
			throw ScoreParamsException ();
		data_->name = const_cast<char*>(matrix.c_str());
		data_->reward = reward;
		data_->penalty = penalty;
		if(Blast_ScoreBlkMatrixFill (data_, 0) != 0)
			throw ScoreParamsException ();
		data_->name = 0;


	}

	// Blastscoreblk(const string &matrix, int gap_open, int gap_extend, int reward, int penalty, const DNA&):
	// data_ (BlastScoreBlkNew(blast_seq_code<DNA>(), 1))
	// {
	// 	data_->name = 0;
	// 	data_->reward = reward;
	// 	data_->penalty = penalty;
	// }


	~Blastscoreblk()
	{ BlastScoreBlkFree(data_); }

/** Supported substitution and gap costs with corresponding quality values
 * for nucleotide sequence comparisons.
 * NB: the values 0 and 0 for the gap costs are treated as the defaults used for
 * the greedy gapped extension, i.e. 
 * gap opening = 0, 
 * gap extension = 1/2 match - mismatch.
 * 
 * The fields are:
 * 
 * 1. Gap opening cost,
 * 2. Gap extension cost,
 * 3. Lambda,
 * 4. K,
 * 5. H,
 * 6. Alpha,
 * 7. Beta,
 * 8. Theta
 */

// TODO: add gumbel parameters for nucleotide cases

/** Karlin-Altschul parameter values for substitution scores 1 and -5. */
// static const array_of_8 blastn_values_1_5[] = {
//     { 0, 0, 1.39, 0.747, 1.38, 1.00,  0, 100 },
//     { 3, 3, 1.39, 0.747, 1.38, 1.00,  0, 100 }
// };


/** Karlin-Altschul parameter values for substitution scores 1 and -4. */
// static const array_of_8 blastn_values_1_4[] = {
//     { 0, 0, 1.383, 0.738, 1.36, 1.02,  0, 100 },
//     { 1, 2,  1.36,  0.67,  1.2,  1.1,  0,  98 }, 
//     { 0, 2,  1.26,  0.43, 0.90,  1.4, -1,  91 },
//     { 2, 1,  1.35,  0.61,  1.1,  1.2, -1,  98 },
//     { 1, 1,  1.22,  0.35, 0.72,  1.7, -3,  88 }
// };

// static const array_of_8 blastn_values_1_1[] = {
//     { 3,  2, 1.09,  0.31, 0.55, 2.0,  -2, 99 },
//     { 2,  2, 1.07,  0.27, 0.49, 2.2,  -3, 97 }, 
//     { 1,  2, 1.02,  0.21, 0.36, 2.8,  -6, 92 }, 
//     { 0,  2, 0.80, 0.064, 0.17, 4.8, -16, 72 },
//     { 4,  1, 1.08,  0.28, 0.54, 2.0,  -2, 98 }, 
//     { 3,  1, 1.06,  0.25, 0.46, 2.3,  -4, 96 }, 
//     { 2,  1, 0.99,  0.17, 0.30, 3.3, -10, 90 }
// };

	template<typename _val>
	int score(_val x, _val y) const
	{ 
		return data_->matrix->data[(long)blast_alphabet<_val>()[(long)AlphabetAttributes<_val>::ALPHABET[x]]][(long)blast_alphabet<_val>()[(long)AlphabetAttributes<_val>::ALPHABET[y]]];
	}
	int score(DNA x, DNA y) const
	{ 
		if (x == y)
		{
			return 4;
		}
		return -5;
		//return getNuclMatchScore((char)AlphabetAttributes<DNA>::ALPHABET[x],(char)AlphabetAttributes<DNA>::ALPHABET[y]);
	}
        // const int blast_lambbda = 1.28;
        // const double blast_k = 0.46;
	//lamda = 0.267

	double lambda() const
	{
		double lamda ;
		if (sequence_type() == amino_acid)
		{
			lamda = data_->kbp_gap_std[0]->Lambda;
		}else if (sequence_type() == nucleotide)
		{
			lamda = 1.09;
		}
		else
		{
			throw ScoreParamsException ();
		}
		return lamda;
	}


	double k() const
	{ 
		double k ;
		if (sequence_type() == amino_acid)
		{
			k = data_->kbp_gap_std[0]->K;
		}else if (sequence_type() == nucleotide)
		{
			k = 0.31;
		}
		else
		{
			throw ScoreParamsException ();
		}
		return k;
	}


	double ln_k() const
	{ 
		double lnk ;
		if (sequence_type() == amino_acid)
		{
			lnk = data_->kbp_gap_std[0]->logK;
		}else if (sequence_type() == nucleotide)
		{
			lnk = log(0.31);
		}
		else
		{
			throw ScoreParamsException ();
		}
		return lnk;
	}


	int low_score() const
	{ 
		int lowscore;
		if (sequence_type() == amino_acid)
		{
			lowscore = data_->loscore;
		
		}else if (sequence_type() == nucleotide)
		{
			lowscore = -4;
		}
		else
		{
			throw ScoreParamsException ();
		}

		// cout<<"lowscore = "<<data_->loscore<<", hisscore = "<<data_->hiscore<<endl;
		return lowscore;
	}

	// void initNuclParas()
	// {
	// 	lamda = 0.267;
	// 	kvalue = 0.041;
	// 	lnk = -3.19;
	// 	lowscore = -4;
	// }

	// void initProteinParas()
	// {
	// 	lamda = data_->kbp_gap_std[0]->Lambda;
	// 	kvalue = data_->kbp_gap_std[0]->K; 
	// 	lnk = data_->kbp_gap_std[0]->logK; 
	// 	lowscore = data_->loscore; 
	// }


private:
	// int kvalue;
	// int lowscore;
	// int lnk;
	// int lamda;
	BlastScoreBlk *data_;

};

const double LN_2 = 0.69314718055994530941723212145818;

class ScoreMatrix
{
	public:
	template<typename _val>
	ScoreMatrix(const string &matrix, int gap_open, int gap_extend, int reward, int penalty, const _val&):
	sb_ (matrix, gap_open, gap_extend, reward, penalty, _val ()),
		bias_ ((char)(-sb_.low_score())),
		name_ (matrix),
		matrix8_ (_val(), sb_),
		matrix8u_ (_val(), sb_, bias_),
		matrix16_ (_val(), sb_)
		
	{ 
		
	}

	// ScoreMatrix(const string &matrix, int gap_open, int gap_extend, int reward, int penalty, const DNA&):
	// 	bias_ ((char)(-sb_.low_score())),
	// 	name_ (matrix),
	// 	matrix8_ (DNA(), sb_),
	// 	matrix8u_ (DNA(), sb_, bias_),
	// 	matrix16_ (DNA(), sb_),
	// 	sb_ (matrix, gap_open, gap_extend, reward, penalty, DNA ())
	// { 

	// }


	template<typename _val>
	void print() const
	{
		// cout << "Scoring matrix = " << name_ << endl;
		// cout << "Lambda = " << sb_.lambda() << endl;
		// cout << "K = " << sb_.k() << endl;
		// const unsigned n = AlphabetAttributes<_val>::ALPHABET_SIZE;
		// for(unsigned i=0;i<n;++i) {
		// 	for(unsigned j=0;j<n;++j)
		// 		printf("%3i", (int)matrix8_.data[i*32+j]);
		// 	printf("\n");
		// }
	}

	static const ScoreMatrix& get()
	{ return *instance; }

	const int8_t* matrix8() const
	{ return matrix8_.data; }

	const uint8_t* matrix8u() const
	{ return matrix8u_.data; }

	const int16_t* matrix16() const
	{ return matrix16_.data; }

	template<typename _val>
	int letter_score(_val a, _val b) const
	{ return matrix8_.data[(int(a) << 5) + int(b)]; }

	template<typename _val>
	uint8_t biased_score(_val a, _val b) const
	{ return matrix8u_.data[(int(a) << 5) + int(b)]; }

	char bias() const
	{ return bias_; }
	
	double bitscore(int raw_score) const
	{ 
		double i = (sb_.lambda() * raw_score - sb_.ln_k() )/ LN_2; 
		//cout<<"bitscore"<<endl;
		// cout<<"bit score = "<<i<<", raw score = "<<raw_score<<", lambda = "<<sb_.lambda()<<", lnk = "<<sb_.ln_k()<<endl;
		return i; 
		//cout<<"bitscore"<<((sb_.lambda() * raw_score - sb_.ln_k() )/ LN_2)<<endl;
	}


	double rawscore(double bitscore, double) const
	{ 
		return (bitscore*LN_2 + sb_.ln_k()) / sb_.lambda(); }


	int rawscore(double bitscore) const
	{ 
		int i =  (int)ceil(rawscore(bitscore, double ())); 
		return i; 
	}

        // double evalue(double bit_score, int db_letters, unsigned query_len) const
        // {
        //     double evalue = log10(query_len) + log10(db_letters) - bit_score * log10(2);
        //     return floor(evalue);
        // }
        // double bitscore(double score) const
        // {
        //     double bit_score = (blast_lambbda * score - log(blast_k)) / log(2);
        //     return bit_score;
        // }

	double evalue(int raw_score, size_t db_letters, unsigned query_len) const
	{ 
		// double evalue = log10(query_len) + log10(db_letters) - raw_score * log10(2);
		double i = static_cast<double>(db_letters) * query_len * pow(2,-bitscore(raw_score)); 
		// cout<<"i = "<<i<<", evale = "<<floor(evalue)<<endl;
		// cout<<"evalue = "<<i<<", raw score = "<<raw_score<<", db letter = "<<db_letters<<", query len = "<<query_len<<endl;
		return i;
	}

	double bitscore(double evalue, size_t db_letters, unsigned query_len) const
	{ 
		return -log(evalue/db_letters/query_len)/log(2); 
	}

	double k() const
	{ return sb_.k(); }


	double lambda() const
	{ return sb_.lambda(); }

	static auto_ptr<ScoreMatrix> instance;

    // const int blast_lambbda = 1.28;
    // const double blast_k = 0.46;
private:

	template<typename _t>
	struct Scores
	{
		template<typename _val>
		Scores(const _val&, const Blastscoreblk &sb, char bias = 0)
		{
			const unsigned n = AlphabetAttributes<_val>::ALPHABET_SIZE;
			for(unsigned i=0;i<32;++i)
				for(unsigned j=0;j<32;++j)
					data[i*32+j] = i < n && j < n ? (_t)(sb.score((_val)i, (_val)j) + (int)bias) : std::numeric_limits<_t>::min();
		}
		
		// template<typename _val>
		Scores(const DNA&, const Blastscoreblk &sb, char bias = 0)
		{
			const unsigned n = AlphabetAttributes<DNA>::ALPHABET_SIZE;
			for(unsigned i=0;i<32;++i)
				for(unsigned j=0;j<32;++j)
					data[i*32+j] = i < n && j < n ? (_t)(sb.score((DNA)i, (DNA)j) + (int)bias) : std::numeric_limits<_t>::min();
		}
		//Tells the compiler to allocate the variable x at a 16-byte aligned memory address instead of the default 4-byte alignment.
		_t data[32*32] __attribute__ ((aligned (16)));
	};

	const Blastscoreblk sb_;
	const char bias_;
	const string name_;
	const Scores<int8_t> matrix8_;
	const Scores<uint8_t> matrix8u_;
	const Scores<int16_t> matrix16_;

};

auto_ptr<ScoreMatrix> ScoreMatrix::instance;



#endif // __NUCLSCOREMATRIX_H__