
#ifndef SEED_H_
#define SEED_H_

#include "VATParameter.h"

typedef uint64_t seed;

unsigned seed_partition(seed s)
{
	return s & (VATParameter::seedp-1);
}

unsigned seed_partition_offset(seed s)
{
	return s >> VATParameter::seedp_bits;
}

#endif /* SEED_H_ */
