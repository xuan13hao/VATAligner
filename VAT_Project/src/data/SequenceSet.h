
#ifndef SEQUENCE_SET_H_
#define SEQUENCE_SET_H_

#include <iostream>
#include <string>
#include "../basic/sequence.h"
#include "AlphabetSet.h"

using std::cout;
using std::endl;
using std::pair;

template<typename _val>
class SequenceSet : public AlphabetSet<_val>
{
	public:
	SequenceSet()
	{ }

	SequenceSet(Input_stream &file):
		AlphabetSet<_val> (file)
	{ }

	void print_stats() const
	{ 
		// cout << "Sequences = " << this->get_length() << ", letters = " << this->letters() << endl; 
		// cout << "Sequences Loading"<< endl; 
	}

	pair<size_t,size_t> len_bounds(size_t min_len) const
	{
		const size_t l (this->get_length());
		size_t max = 0, min = std::numeric_limits<size_t>::max();
		for(size_t i=0;i<l;++i) {
			max = std::max(this->length(i), max);
			min = this->length(i) >= min_len ? std::min(this->length(i), min) : min;
		}
		return pair<size_t,size_t> (min, max);
	}

	sequence<const _val> window_infix(size_t offset, unsigned &left) const
	{
		const _val* begin (this->data(offset));
		unsigned n (0);
		while(*begin != AlphabetSet<_val>::PADDING_CHAR && n <= VATParameters::window) {
			--begin;
			++n;
		}
		++begin;
		left = VATParameters::window + 1 - n;
		const _val* end (this->data(offset));
		n = 0;
		while(*end != AlphabetSet<_val>::PADDING_CHAR && n < VATParameters::window) {
			++end;
			++n;
		}
		return sequence<const _val> (begin, end - begin);
	}

	sequence<const _val> fixed_window_infix(size_t offset) const
	{
		const _val* begin (this->data(offset));
		unsigned n (0);
		while(*begin != AlphabetSet<_val>::PADDING_CHAR && n <= VATParameters::window) {
			--begin;
			++n;
		}
		++begin;
		const _val* s (this->data(offset - VATParameters::window));
		return sequence<const _val> (s, 2*VATParameters::window, begin - s);
	}

	vector<size_t> partition() const
	{
		vector<size_t> v;
		const size_t l = (this->letters()+VATConsts::seqp-1) / VATConsts::seqp;
		v.push_back(0);
		for(unsigned i=0;i<this->get_length();) {
			size_t n = 0;
			while(i<this->get_length() && n < l)
				n += this->length(i++);
			v.push_back(i);
		}
		for(unsigned i=v.size();i<VATConsts::seqp+1;++i)
			v.push_back(this->get_length());
		
		// for(int i = 0; i <v.size();i++)
		// {
		// 	cout<<"partition = "<<v[i]<<endl;
		// }
		return v;
	}

	size_t reverse_translated_len(size_t i) const
	{
		const size_t j (i - i%6);
		const size_t l (this->length(j));
		if(this->length(j+2) == l)
			return l*3 + 2;
		else if(this->length(j+1) == l)
			return l*3 + 1;
		else
			return l*3;
	}

	virtual ~SequenceSet()
	{ }

};

#endif /* SEQUENCE_SET_H_ */
