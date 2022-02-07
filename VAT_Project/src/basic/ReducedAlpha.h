
#ifndef REDUCTION_H_
#define REDUCTION_H_

using std::string;
using std::vector;

#include "value.h"

typedef enum { dna = 0, protein=1, rna =2 } ReducedAlphaType;

template<typename _val>
class ReducedAlpha
{
	public:
    ReducedAlpha(int i)
	{
		memset(map_, 0, sizeof(map_));
		memset(map8_, 0, sizeof(map8_));
		switch (i)
		{
			case 0:
				initDNA();
				break;
			case 1:
				initProtein();
				break;
			case 2:
				initRNA();
				break;

		}
		
	}
//KREDQN C G H M F Y ILV W P STA
	void initProtein()
	{
		map_[(long)AlphabetFeature<_val>::bioa('K')] = 0;
		map_[(long)AlphabetFeature<_val>::bioa('R')] = 0;
		map_[(long)AlphabetFeature<_val>::bioa('E')] = 0;
		map_[(long)AlphabetFeature<_val>::bioa('D')] = 0;
		map_[(long)AlphabetFeature<_val>::bioa('Q')] = 0;
		map_[(long)AlphabetFeature<_val>::bioa('N')] = 0;

		map_[(long)AlphabetFeature<_val>::bioa('C')] = 1;
		map_[(long)AlphabetFeature<_val>::bioa('G')] = 2;
		map_[(long)AlphabetFeature<_val>::bioa('H')] = 3;
		map_[(long)AlphabetFeature<_val>::bioa('M')] = 4;
		map_[(long)AlphabetFeature<_val>::bioa('F')] = 5;
		map_[(long)AlphabetFeature<_val>::bioa('Y')] = 6;

		map_[(long)AlphabetFeature<_val>::bioa('I')] = 7;
		map_[(long)AlphabetFeature<_val>::bioa('L')] = 7;
		map_[(long)AlphabetFeature<_val>::bioa('V')] = 7;

		map_[(long)AlphabetFeature<_val>::bioa('W')] = 8;
		map_[(long)AlphabetFeature<_val>::bioa('P')] = 9;

		map_[(long)AlphabetFeature<_val>::bioa('S')] = 10;
		map_[(long)AlphabetFeature<_val>::bioa('T')] = 10;
		map_[(long)AlphabetFeature<_val>::bioa('A')] = 10;

		map8_[(long) AlphabetFeature<_val>::bioa('K')] =  0;
		map8_[(long) AlphabetFeature<_val>::bioa('R')] =  0;
		map8_[(long) AlphabetFeature<_val>::bioa('E')] =  0;
		map8_[(long) AlphabetFeature<_val>::bioa('D')] =  0;
		map8_[(long) AlphabetFeature<_val>::bioa('Q')] =  0;
		map8_[(long) AlphabetFeature<_val>::bioa('N')] =  0;

		map8_[(long) AlphabetFeature<_val>::bioa('C')] =  1;
		map8_[(long) AlphabetFeature<_val>::bioa('G')] =  2;
		map8_[(long) AlphabetFeature<_val>::bioa('H')] =  3;
		map8_[(long) AlphabetFeature<_val>::bioa('M')] =  4;
		map8_[(long) AlphabetFeature<_val>::bioa('F')] =  5;
		map8_[(long) AlphabetFeature<_val>::bioa('Y')] =  6;

		map8_[(long) AlphabetFeature<_val>::bioa('I')] =  7;
		map8_[(long) AlphabetFeature<_val>::bioa('L')] =  7;
		map8_[(long) AlphabetFeature<_val>::bioa('Y')] =  7;

		map8_[(long) AlphabetFeature<_val>::bioa('W')] =  8;
		map8_[(long) AlphabetFeature<_val>::bioa('P')] =  9;

		map8_[(long) AlphabetFeature<_val>::bioa('S')] =  10;
		map8_[(long) AlphabetFeature<_val>::bioa('T')] =  10;
		map8_[(long) AlphabetFeature<_val>::bioa('A')] =  10;

	}

	void initDNA()
	{
		map_[(long)AlphabetFeature<_val>::bioa('A')] = 0;
		map8_[(long) AlphabetFeature<_val>::bioa('A')] =  0;
		map_[(long)AlphabetFeature<_val>::bioa('C')] = 1;
		map8_[(long) AlphabetFeature<_val>::bioa('C')] =  1;
		map_[(long)AlphabetFeature<_val>::bioa('G')] = 2;
		map8_[(long) AlphabetFeature<_val>::bioa('G')] =  2;
		map_[(long)AlphabetFeature<_val>::bioa('T')] = 3;
		map8_[(long) AlphabetFeature<_val>::bioa('T')] =  3;
	}

	void initRNA()
	{
		map_[(long)AlphabetFeature<_val>::bioa('A')] = 0;
		map8_[(long) AlphabetFeature<_val>::bioa('A')] =  0;
		map_[(long)AlphabetFeature<_val>::bioa('C')] = 1;
		map8_[(long) AlphabetFeature<_val>::bioa('C')] =  1;
		map_[(long)AlphabetFeature<_val>::bioa('G')] = 2;
		map8_[(long) AlphabetFeature<_val>::bioa('G')] =  2;
		map_[(long)AlphabetFeature<_val>::bioa('U')] = 3;
		map8_[(long) AlphabetFeature<_val>::bioa('U')] =  3;
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

template<> const ReducedAlpha<Protein> ReducedAlpha<Protein>::reduction (ReducedAlphaType::protein);
template<> const ReducedAlpha<DNA> ReducedAlpha<DNA>::reduction (ReducedAlphaType::dna);
template<> const ReducedAlpha<RNA> ReducedAlpha<RNA>::reduction (ReducedAlphaType::rna);


#ifdef EXTRA
#include "../../../extra/reduction.h"
#endif

#endif /* REDUCTION_H_ */
