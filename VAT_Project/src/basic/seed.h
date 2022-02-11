
#ifndef SEED_H_
#define SEED_H_

#include "const.h"

typedef uint64_t seed;

unsigned seed_partition(seed s)
{
	return s & (Const::seedp-1);
}

unsigned seed_partition_offset(seed s)
{
	return s >> Const::seedp_bits;
}

#endif /* SEED_H_ */
