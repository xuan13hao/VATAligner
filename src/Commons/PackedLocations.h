

#ifndef PACKED_LOC_H_
#define PACKED_LOC_H_

class packed_uint40_t
{
	public:
	uint8_t		high;
	uint32_t	low;
	packed_uint40_t():
		high (),
		low ()
	{ }
	packed_uint40_t(uint64_t v):
		high (v>>32),
		low (v&0xffffffffu)
	{ }
	operator const uint64_t() const
	{ return (uint64_t(high) << 32) | low; }
	bool operator<(const packed_uint40_t &rhs) const
	{ return high < rhs.high || (high == rhs.high && low < rhs.low); }
	friend uint64_t operator-(const packed_uint40_t &x, const packed_uint40_t &y)
	{ return (const uint64_t)(x) - (const uint64_t)(y); }
} __attribute__((packed));

template<class _loc>
struct packed_sequence_location
{
	typedef _loc type;
};

template<>
struct packed_sequence_location<uint64_t>
{
	typedef packed_uint40_t type;
};

#endif /* PACKED_LOC_H_ */
