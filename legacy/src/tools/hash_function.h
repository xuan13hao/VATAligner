

#ifndef HASH_FUNCTION_H_
#define HASH_FUNCTION_H_

class murmur_hash {
	public:
	
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

class hash64 {
	public:
	uint64_t operator()(uint64_t key, uint64_t mask)
	{
		key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
		key = key ^ key >> 24;
		key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
		key = key ^ key >> 14;
		key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
		key = key ^ key >> 28;
		key = (key + (key << 31)) & mask;
		return key;
	}
};
#endif /* HASH_FUNCTION_H_ */
