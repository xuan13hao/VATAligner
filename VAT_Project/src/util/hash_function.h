

#ifndef HASH_FUNCTION_H_
#define HASH_FUNCTION_H_

struct murmur_hash {
	uint64_t operator()(uint64_t h) const
	{
		h ^= h >> 33;
		h *= 0xff51afd7ed558ccdLL;
		h ^= h >> 33;
		h *= 0xc4ceb9fe1a85ec53LL;
		h ^= h >> 33;
		return h;
	}
};

#endif /* HASH_FUNCTION_H_ */
