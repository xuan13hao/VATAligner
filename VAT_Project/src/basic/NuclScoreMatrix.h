#ifndef __NUCLSCOREMATRIX_H__
#define __NUCLSCOREMATRIX_H__


#include "../algo/blast/core/blast_stat.h"
#include "../algo/blast/core/blast_encoding.h"
#include "ProteinProfile.h"
#include "Matrixs.h"
#include "AlphabetType.h"
using std::string;
using std::cout;
using std::endl;
using std::auto_ptr;

class ScoreParamsException : public std::exception
{
	public:
	virtual const char* what() const throw()
	{ 
		return "Error scoring parameters"; 
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
			throw Score_params_exception ();
		if((data_->kbp_gap_std[0] = Blast_KarlinBlkNew()) == 0
				|| (data_->kbp_std[0] = Blast_KarlinBlkNew()) == 0)
			throw Score_params_exception ();
		if(blast_load_karlin_blk<_val>(data_->kbp_gap_std[0],
				data_->kbp_std[0],
				gap_open,
				gap_extend,
				reward,
				penalty,
				matrix.c_str()) != 0)
			throw Score_params_exception ();
		data_->name = const_cast<char*>(matrix.c_str());
		data_->reward = reward;
		data_->penalty = penalty;
		if(Blast_ScoreBlkMatrixFill (data_, 0) != 0)
			throw Score_params_exception ();
		data_->name = 0;
	}


	Blastscoreblk(const string &matrix, int gap_open, int gap_extend, int reward, int penalty, const DNA&)
	{
		data_->name = 0;
		data_->reward = reward;
		data_->penalty = penalty;
	}



	~Blastscoreblk()
	{ BlastScoreBlkFree(data_); }


	template<typename _val>
	int score(_val x, _val y) const
	{ 
		return data_->matrix->data[(long)blast_alphabet<_val>()[(long)AlphabetAttributes<_val>::ALPHABET[x]]][(long)blast_alphabet<_val>()[(long)AlphabetAttributes<_val>::ALPHABET[y]]];
	}
	int score(DNA x, DNA y) const
	{ 
		if (x == y)
		{
			return 1;
		}
		return 0;
		//return getNuclMatchScore((char)AlphabetAttributes<DNA>::ALPHABET[x],(char)AlphabetAttributes<DNA>::ALPHABET[y]);
	}
        // const int blast_lambbda = 1.28;
        // const double blast_k = 0.46;
	//lamda = 0.267
	template<typename _val>
	double lambda(_val&) const
	{ 
		double lamda = 0.267;
		return lamda;
	}

	double lambda(DNA()) const
	{ 
		double lamda = 0.267;
		return lamda;
	}
	//k = 0.041
	template<typename _val>
	double k(_val&) const
	{ 
		double k = 0.041;
		return k;
	}
	double k(DNA()) const
	{ 
		double k = 0.041;
		return k;
	}
	template<typename _val>
	double ln_k(_val&) const
	{ 
		double lnk = -3.19;
		return lnk;
	}
	double ln_k(DNA()) const
	{ 
		double lnk = -3.19;
		return lnk;
	}

	template<typename _val>
	int low_score(_val&) const
	{ 
		// cout<<"low_score = "<<data_->loscore<<endl;
		int lowscore = -4;
		return lowscore;
	}
	int low_score(DNA()) const
	{ 
		// cout<<"low_score = "<<data_->loscore<<endl;
		int lowscore = -4;
		return lowscore;
	}
private:

	BlastScoreBlk *data_;

};

const double LN_2 = 0.69314718055994530941723212145818;

class NuclScoreMatrix
{
	public:
	template<typename _val>
	NuclScoreMatrix(const string &matrix, int gap_open, int gap_extend, int reward, int penalty, const _val&):
		sb_ (matrix, gap_open, gap_extend, reward, penalty, _val ()),
		bias_ ((char)(-sb_.low_score())),
		name_ (matrix),
		matrix8_ (_val(), sb_),
		matrix8u_ (_val(), sb_, bias_),
		matrix16_ (_val(), sb_)
	{ }
/*
	score_matrix(const string &matrix, int gap_open, int gap_extend, int reward, int penalty, const DNA&):
		sb_ (matrix, gap_open, gap_extend, reward, penalty, DNA ()),
		bias_ ((char)(-sb_.low_score())),
		name_ (matrix),
		matrix8_ (DNA(), sb_),
		matrix8u_ (DNA(), sb_, bias_),
		matrix16_ (DNA(), sb_)
	{ 
		
	}

*/
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

	static const NuclScoreMatrix& get()
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
	template<typename _val>
	double bitscore(int raw_score) const
	{ return ( sb_.lambda(_val&) * raw_score - sb_.ln_k(_val&)) / LN_2; }

	template<typename _val>
	double rawscore(double bitscore, double) const
	{ return (bitscore*LN_2 + sb_.ln_k(_val&)) / sb_.lambda(_val&); }
	template<typename _val>
	int rawscore(double bitscore,_val&) const
	{ 
		int i =  (int)ceil(rawscore(bitscore, double ())); 
		return i; 
	}
	template<typename _val>
	double evalue(int raw_score, size_t db_letters, unsigned query_len) const
	{ return static_cast<double>(db_letters) * query_len * pow(2,-bitscore(raw_score)); }

	double bitscore(double evalue, size_t db_letters, unsigned query_len) const
	{ return -log(evalue/db_letters/query_len)/log(2); }

	template<typename _val>
	double k() const
	{ return sb_.k(_val&); }

	template<typename _val>
	double lambda() const
	{ return sb_.lambda(_val&); }

	static auto_ptr<NuclScoreMatrix> instance;

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
		
		_t data[32*32] __attribute__ ((aligned (16)));
	};

	const Blastscoreblk sb_;
	const char bias_;
	const string name_;
	const Scores<int8_t> matrix8_;
	const Scores<uint8_t> matrix8u_;
	const Scores<int16_t> matrix16_;

};

auto_ptr<NuclScoreMatrix> NuclScoreMatrix::instance;



#endif // __NUCLSCOREMATRIX_H__