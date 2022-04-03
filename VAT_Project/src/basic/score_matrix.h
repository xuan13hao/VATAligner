
#ifndef SCORE_MATRIX_H_
#define SCORE_MATRIX_H_

#include "../algo/blast/core/blast_stat.h"
#include "../algo/blast/core/blast_encoding.h"
#include "score_traits.h"
#include "Matrixs.h"
using std::string;
using std::cout;
using std::endl;
using std::auto_ptr;

struct Score_params_exception : public std::exception
{
	virtual const char* what() const throw()
	{ return "Invalid scoring parameters"; }
};

struct Blast_score_blk
{

	template<typename _val>
	Blast_score_blk(const string &matrix, int gap_open, int gap_extend, int reward, int penalty, const _val&):
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

	~Blast_score_blk()
	{ BlastScoreBlkFree(data_); }
	// for protein alignment
	template<typename _val>
	int score(_val x, _val y) const
	{ 

		return data_->matrix->data[(long)blast_alphabet<_val>()[(long)AlphabetAttributes<_val>::ALPHABET[x]]][(long)blast_alphabet<_val>()[(long)AlphabetAttributes<_val>::ALPHABET[y]]];
	}
	//template<typename _val>
	
	int score(DNA x, DNA y) const
	{ 
		return getNuclMatchScore((char)AlphabetAttributes<DNA>::ALPHABET[x],(char)AlphabetAttributes<DNA>::ALPHABET[y]);
	}

	//lamda = 0.267
	double lambda() const
	{ 
		// double lamda = 0.267;
		// return lamda;
		// cout<<"lamda = "<<data_->kbp_gap_std[0]->Lambda<<endl;
		return data_->kbp_gap_std[0]->Lambda; 
	}
	//k = 0.041
	double k() const
	{ 
		// double k = 0.041;
		// return k;
		// cout<<"k = "<<data_->kbp_gap_std[0]->K<<endl;
		return data_->kbp_gap_std[0]->K; 
	}
	//lnk = -3.19
	double ln_k() const
	{ 
		// double lnk = -3.19;
		// return lnk;
		// cout<<"ln k = "<<data_->kbp_gap_std[0]->logK<<endl;
		return data_->kbp_gap_std[0]->logK; 
	}
	//low socre = -4
	int low_score() const
	{ 
		// int lowscore = -4;
		// return lowscore;
		//cout<<"low_score = "<<data_->loscore<<endl;
		return data_->loscore; 
	}

private:

	BlastScoreBlk *data_;

};

const double LN_2 = 0.69314718055994530941723212145818;

struct score_matrix
{

	template<typename _val>
	score_matrix(const string &matrix, int gap_open, int gap_extend, int reward, int penalty, const _val&):
		sb_ (matrix, gap_open, gap_extend, reward, penalty, _val ()),
		bias_ ((char)(-sb_.low_score())),
		name_ (matrix),
		matrix8_ (_val(), sb_),
		matrix8u_ (_val(), sb_, bias_),
		matrix16_ (_val(), sb_)
	{ }

	score_matrix(const string &matrix, int gap_open, int gap_extend, int reward, int penalty, const DNA&):
		sb_ (matrix, gap_open, gap_extend, reward, penalty, DNA ()),
		bias_ ((char)(-sb_.low_score())),
		name_ (matrix),
		matrix8_ (DNA(), sb_),
		matrix8u_ (DNA(), sb_, bias_),
		matrix16_ (DNA(), sb_)
	{ 
		
	}


	template<typename _val>
	void print() const
	{
		cout << "Scoring matrix = " << name_ << endl;
		// cout << "Lambda = " << sb_.lambda() << endl;
		// cout << "K = " << sb_.k() << endl;
		// const unsigned n = AlphabetAttributes<_val>::ALPHABET_SIZE;
		// for(unsigned i=0;i<n;++i) {
		// 	for(unsigned j=0;j<n;++j)
		// 		printf("%3i", (int)matrix8_.data[i*32+j]);
		// 	printf("\n");
		// }
	}

	static const score_matrix& get()
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
	{ return ( sb_.lambda() * raw_score - sb_.ln_k()) / LN_2; }

	double rawscore(double bitscore, double) const
	{ return (bitscore*LN_2 + sb_.ln_k()) / sb_.lambda(); }

	int rawscore(double bitscore) const
	{ 
		// return (int)ceil(rawscore(bitscore, double ())); 
		return 10; 
	}

	double evalue(int raw_score, size_t db_letters, unsigned query_len) const
	{ return static_cast<double>(db_letters) * query_len * pow(2,-bitscore(raw_score)); }

	double bitscore(double evalue, size_t db_letters, unsigned query_len) const
	{ return -log(evalue/db_letters/query_len)/log(2); }

	double k() const
	{ return sb_.k(); }

	double lambda() const
	{ return sb_.lambda(); }

	static auto_ptr<score_matrix> instance;

private:

	template<typename _t>
	struct Scores
	{
		template<typename _val>
		Scores(const _val&, const Blast_score_blk &sb, char bias = 0)
		{
			const unsigned n = AlphabetAttributes<_val>::ALPHABET_SIZE;
			for(unsigned i=0;i<32;++i)
				for(unsigned j=0;j<32;++j)
					data[i*32+j] = i < n && j < n ? (_t)(sb.score((_val)i, (_val)j) + (int)bias) : std::numeric_limits<_t>::min();
		}
		// template<typename _val>
		Scores(const DNA&, const Blast_score_blk &sb, char bias = 0)
		{
			const unsigned n = AlphabetAttributes<DNA>::ALPHABET_SIZE;
			for(unsigned i=0;i<32;++i)
				for(unsigned j=0;j<32;++j)
					data[i*32+j] = i < n && j < n ? (_t)(sb.score((DNA)i, (DNA)j) + (int)bias) : std::numeric_limits<_t>::min();
		}
		
		_t data[32*32] __attribute__ ((aligned (16)));
	};

	const Blast_score_blk sb_;
	const char bias_;
	const string name_;
	const Scores<int8_t> matrix8_;
	const Scores<uint8_t> matrix8u_;
	const Scores<int16_t> matrix16_;

};

auto_ptr<score_matrix> score_matrix::instance;

#endif /* SCORE_MATRIX_H_ */
