
#ifndef REDUCTION_H_
#define REDUCTION_H_

using std::string;
using std::vector;

#include "value.h"

template<typename _val>
class ReducedAlpha
{
	public:
    ReducedAlpha(const char *definition_string)
	{
		memset(map_, 0, sizeof(map_));
		memset(map8_, 0, sizeof(map8_));
		const vector<string> tokens (tokenize(definition_string, " "));
		size_ = tokens.size();
        for(unsigned i=0;i<size_;++i)
        	for(unsigned j=0;j<tokens[i].length();++j) {
        		const char ch = tokens[i][j];
        		map_[(long) AlphabetAttributes<_val>::from_char(ch)] =  i;
                map8_[(long) AlphabetAttributes<_val>::from_char(ch)] =  i;
            }
	}

	unsigned size() const
	{ return size_; }

	unsigned operator()(_val a) const
	{ return map_[(long)a]; }

	const char* map8() const
	{ return map8_; }

	static const ReducedAlpha reduction;

private:

	unsigned map_[256];
	char map8_[256];
	unsigned size_;

};
/*
//hsdm.4
char hsdm4s[500][20] = {"LIVFMYW", "C", "DNTSKEQRAGP", "H"};
DIAMOND 11
KREDQN C G H M F Y ILV W P STA"
*/
// template<> const ReducedAlpha<Protein> ReducedAlpha<Protein>::reduction ("KREDQN C G H M F Y ILV W P STA");
template<> const ReducedAlpha<Protein> ReducedAlpha<Protein>::reduction ("LIVFMYW C DNTSKEQRAGP H");
template<> const ReducedAlpha<DNA> ReducedAlpha<DNA>::reduction ("A C G T");


#endif /* REDUCTION_H_ */
