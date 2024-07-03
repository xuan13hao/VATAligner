
#ifndef SEED_H_
#define SEED_H_

#include "VATConsts.h"

typedef uint64_t seed;

unsigned seed_partition(seed s)
{
	return s & (VATConsts::seedp-1);
}

unsigned seed_partition_offset(seed s)
{
	// return s & 1024;
	return s >> VATConsts::seedp_bits;
}

#endif /* SEED_H_ */
