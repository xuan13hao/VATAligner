

#ifndef SHAPE_H_
#define SHAPE_H_

#include "VATConsts.h"
#include "value.h"
#include "SeedPartition.h"
#include "../tools/hash_function.h"
#include "NuclScoreMatrix.h"
#include "ReducedAlpha.h"

// const double background_freq[] = {-1.188861,
// 		-4.343446,
// 		-2.648093,
// 		-3.806941,
// 		-3.742636,
// 		-3.221182,
// 		-3.498273,
// 		-1.498637,
// 		-4.339607,
// 		-3.027002,
// 		-1.557546 };


// const double background_freq[] = {-1.38,-1.38,-1.38,-1.38};

const double background_freq[] = {-3,-3,-3,-3};
template<typename _val>
bool use_seed_freq()
{ return false; }

template<>
bool use_seed_freq<Protein>()
{ return true; }

template<>
bool use_seed_freq<DNA>()
{ return true; }


class Shape
{
	public:
	Shape():
		length_ (0),
		weight_ (0),
		d_ (0),
		mask_ (0),
		rev_mask_ (0),
		id_ (0)
	{ memset(positions_, 0, sizeof(uint32_t)*VATConsts::max_seed_weight); }

	Shape(const char *code, unsigned id):
		weight_ (0),
		mask_ (0),
		rev_mask_ (0),
		id_ (id)
	{
		assert(id < VATConsts::max_shapes);
		assert(strlen(code) <= 32);
		memset(positions_, 0, sizeof(uint32_t)*VATConsts::max_seed_weight);
		unsigned i (0);
		for(;i<strlen(code);++i) 
		{
			rev_mask_ <<= 1;
			if(code[i] == '1') 
			{
				assert(weight_ < VATConsts::max_seed_weight);
				positions_[weight_] = i;
				++weight_;
				mask_ |= 1 << i;
				rev_mask_ |= 1;
			}
		}
		length_ = i;
		d_ = positions_[weight_/2-1];
	}
/**
* shape with minimizer
*/



	
	//template<>
	inline bool set_seed(uint64_t &s, const Protein *seq) const
	{
		s = 0;
		double f = 0;
		for(unsigned i=0;i<weight_;++i) 
		{
			Protein l = seq[positions_[i]];
			if(l == AlphabetAttributes<Protein>::MASK_CHAR || l == AlphabetSet<Protein>::PADDING_CHAR)
				return false;
			l = mask_critical(l);
			unsigned r = ReducedAlpha<Protein>::reduction(l);
			f += background_freq[r];
			s *= ReducedAlpha<Protein>::reduction.size();
			s += uint64_t(r);
		}
		s = murmur_hash()(s);
		if(use_seed_freq<Protein>() && f > VATParameters::max_seed_freq) return false;

		return true;
	}

	inline bool set_seed(uint64_t &s, const DNA *seq) const
	{
		s = 0;
		double f = 0;
		for(unsigned i=0;i<weight_;++i) 
		{
			DNA l = seq[positions_[i]];
			if(l == AlphabetAttributes<DNA>::MASK_CHAR || l == AlphabetSet<DNA>::PADDING_CHAR)
				return false;
			l = mask_critical(l);
			unsigned r = ReducedAlpha<DNA>::reduction(l);
			f += background_freq[r];
			s *= ReducedAlpha<DNA>::reduction.size();
			s += uint64_t(r);
		}
		s = murmur_hash()(s);
		if(use_seed_freq<DNA>() && f > VATParameters::max_seed_freq) return false;

		return true;
	}

inline uint64_t calculate_minimizer(const DNA *seq, unsigned start, unsigned length)
{
    const unsigned k = length;  // Size of the minimizer
    const uint64_t max_val = std::numeric_limits<uint64_t>::max();

    uint64_t min_seed = max_val;
    for (unsigned i = start; i < start + length - k + 1; ++i)
    {
        uint64_t seed;
        if (set_seed(seed, &seq[i]))
        {
            if (seed < min_seed)
                min_seed = seed;
        }
    }

    return min_seed;
}



	//template<>
	inline bool set_seed_minimizer(uint64_t &s, const DNA *seq) const
	{

		s = 0;
		double f = 0;
		uint64_t minimizer = 0;
		uint64_t mask = (1ULL<<2*weight_) - 1;
		for(unsigned i=0;i<weight_;++i) 
		{
			DNA l = seq[positions_[i]];
			if(l == AlphabetAttributes<DNA>::MASK_CHAR || l == AlphabetSet<DNA>::PADDING_CHAR)
				return false;
			l = mask_critical(l);
			unsigned r = ReducedAlpha<DNA>::reduction(l);
			// cout<<"r = "<<r<<endl;
			minimizer = minimizer + r*pow(4,weight_-i-1);
			// cout<<"k-1 = "<<weight_-i-1<<",r = "<<r<<", 4^k-1 = "<<pow(4,weight_-i-1)<<", minimizer = "<<minimizer<<endl;
			f += background_freq[r];
			// s *= ReducedAlpha<_val>::reduction.size();
			// cout<<AlphabetAttributes<_val>::ALPHABET[r];
			// s += uint64_t(r);
		}
		// cout<<endl;
		s = hash64()(minimizer,mask);
		// s = murmur_hash()(s);
		// cout<<"minimizer = "<<minimizer<<endl;
		if(use_seed_freq<DNA>() && f > VATParameters::max_seed_freq) return false;

		return true;
	}

	template<typename _val>
	inline bool	is_low_freq(const _val *seq) const
	{
		double f = 0;
		for(unsigned i=0;i<weight_;++i) {
			_val l = seq[positions_[i]];
			if(l == AlphabetAttributes<_val>::MASK_CHAR || l == AlphabetSet<_val>::PADDING_CHAR)
				return false;
			l = mask_critical(l);
			unsigned r = ReducedAlpha<_val>::reduction(l);
			f += background_freq[r];
		}
		return !use_seed_freq<_val>() || f <= VATParameters::max_seed_freq;
	}

	template<typename _val>
	inline bool	is_low_freq_rev(const _val *seq) const
	{
		double f = 0;
		for(unsigned i=0;i<weight_;++i) {
			_val l = seq[(int)positions_[i]-(int)length_];
			if(l == AlphabetAttributes<_val>::MASK_CHAR || l == AlphabetSet<_val>::PADDING_CHAR)
				return false;
			l = mask_critical(l);
			unsigned r = ReducedAlpha<_val>::reduction(l);
			f += background_freq[r];
		}
		return !use_seed_freq<_val>() || f <= VATParameters::max_seed_freq;
	}

	uint32_t length_, weight_, positions_[VATConsts::max_seed_weight], d_, mask_, rev_mask_, id_;

};

#endif /* SHAPE_H_ */
