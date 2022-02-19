
#ifndef REDUCTION_H_
#define REDUCTION_H_

using std::string;
using std::vector;

#include "value.h"

template<typename _val>
struct Reduction
{

    Reduction(const char *definition_string)
	{
		memset(map_, 0, sizeof(map_));
		memset(map8_, 0, sizeof(map8_));
		const vector<string> tokens (tokenize(definition_string, " "));
		size_ = tokens.size();
        for(unsigned i=0;i<size_;++i)
        	for(unsigned j=0;j<tokens[i].length();++j) {
        		const char ch = tokens[i][j];
        		map_[(long) Value_traits<_val>::from_char(ch)] =  i;
                map8_[(long) Value_traits<_val>::from_char(ch)] =  i;
            }
	}

	unsigned size() const
	{ return size_; }

	unsigned operator()(_val a) const
	{ return map_[(long)a]; }

	const char* map8() const
	{ return map8_; }

	static const Reduction reduction;

private:

	unsigned map_[256];
	char map8_[256];
	unsigned size_;

};

template<> const Reduction<Protein> Reduction<Protein>::reduction ("KREDQN C G H M F Y ILV W P STA");
template<> const Reduction<DNA> Reduction<DNA>::reduction ("A C G T");


#endif /* REDUCTION_H_ */
