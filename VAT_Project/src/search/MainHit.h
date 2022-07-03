

#ifndef COLLISION_H_
#define COLLISION_H_

template<typename _val>
inline unsigned letter_match(const _val query, const _val subject)
{
	if(query != AlphabetAttributes<_val>::MASK_CHAR && ReducedAlpha<_val>::reduction(query) == ReducedAlpha<_val>::reduction(mask_critical(subject)))
		return 1;
	else
		return 0;
}

inline bool match_shape_mask(const uint64_t mask, const uint64_t shape_mask)
{
	return (mask & shape_mask) == shape_mask;
}

template<typename _val>
inline bool is_lower_chunk(const _val *subject, unsigned sid)
{
	uint64_t seed;
	ShapeConfigures::get().get_shape(sid).set_seed(seed, subject);
	return current_range.lower(seed_partition(seed));
}

template<typename _val>
inline bool is_lower_or_equal_chunk(const _val *subject, unsigned sid)
{
	uint64_t seed;
	ShapeConfigures::get().get_shape(sid).set_seed(seed, subject);
	return current_range.contains(seed_partition(seed));
}

inline bool need_lookup(unsigned sid)
{
	//return sid != 0 || current_query_chunk != 0;
	return sid != 0;
}

template<typename _val>
inline bool is_low_freq(const _val *subject, unsigned sid)
{ return ShapeConfigures::get().get_shape(sid).is_low_freq(subject); }

template <typename _val, typename _pos>
inline bool shape_collision_right(uint64_t mask, uint64_t shape_mask, const _val *subject, unsigned sid)
{
	// bool c = (is_lower_chunk(subject, sid));
	// // bool g = (need_lookup(sid) && (!ReferenceSeqs<_val>::get().get_masking(subject, sid)));
	// cout<<"right"<<", c = "<<c<<endl;

	if(!match_shape_mask(mask, shape_mask)) return false;
	return is_lower_chunk(subject, sid)
			&& is_low_freq(subject, sid)
			&& (!get_critical(*subject) || (need_lookup(sid) && !ReferenceSeqs<_val>::get().get_masking(subject, sid)));
}

template <typename _val, typename _pos>
inline bool shape_collision_left(uint64_t mask, uint64_t shape_mask, const _val *subject, unsigned sid, bool chunked)
{
	// bool c = (is_lower_or_equal_chunk(subject, sid));
	// // bool g = (need_lookup(sid) && (!ReferenceSeqs<_val>::get().get_masking(subject, sid)));
	// cout<<"left "<<", c = "<<c<<endl;
	if(!match_shape_mask(mask, shape_mask))
		return false;
	return (!chunked && (!get_critical(*subject) || (need_lookup(sid)&& !ReferenceSeqs<_val>::get().get_masking(subject, sid)))
	&& is_low_freq(subject, sid)
	);
	// return (!chunked || is_lower_or_equal_chunk(subject, sid))
	// 		&& is_low_freq(subject, sid)
	// 		&& (!get_critical(*subject) || (need_lookup(sid) && !ReferenceSeqs<_val>::get().get_masking(subject, sid)));
}

template <typename _val, typename _pos>
inline bool previous_shape_collision(uint64_t mask, uint64_t shape_mask, const _val *subject, unsigned sid)
{
	if(!match_shape_mask(mask, shape_mask)) return false;
	return is_low_freq(subject, sid)
			&& (!get_critical(*subject) || !ReferenceSeqs<_val>::get().get_masking(subject, sid));
}

template<typename _val, typename _pos>
bool is_primary_hit2(const _val *query,
					const _val *subject,
					const unsigned seed_offset,
					const unsigned sid,
					const unsigned len,
					Statistics &stat)
{
	assert(len > 0 && len <= VATParameters::window*2);
	unsigned mask (0);
	const bool chunked (VATParameters::lowmem > 1);
	unsigned i = len;

	do {
		--i;
		mask <<= 1;
		mask |= letter_match(query[i], subject[i]);
		for(unsigned j=0;j<sid;++j)
			if(previous_shape_collision<_val,_pos>(mask, ShapeConfigures::instance.get_shape(j).mask_, &subject[i], j, stat))
				return false;
		if(chunked && i > seed_offset && shape_collision_right<_val,_pos>(mask, ShapeConfigures::instance.get_shape(sid).mask_, &subject[i], sid, stat))
			return false;
	} while (i > seed_offset);

	if(i == 0)
		return true;

	do {
		--i;
		mask <<= 1;
		mask |= letter_match(query[i], subject[i]);
		for(unsigned j=0;j<sid;++j)
			if(previous_shape_collision<_val,_pos>(mask, ShapeConfigures::instance.get_shape(j).mask_, &subject[i], j, stat))
				return false;
		if(shape_collision_left<_val,_pos>(mask, ShapeConfigures::instance.get_shape(sid).mask_, &subject[i], sid, chunked, stat))
			return false;
	} while (i > 0);

	return true;
}

template<typename _val, typename _pos>
bool is_primary_hit(const _val *query,
					const _val *subject,
					const unsigned seed_offset,
					const unsigned sid,
					const unsigned len)
{
	assert(len > 0 && len <= VATParameters::window*2);
	const bool chunked (VATParameters::lowmem > 1);

	uint64_t mask = reduced_match32(query, subject, len);
	unsigned i = 0;
	uint64_t current_mask = ShapeConfigures::instance.get_shape(sid).mask_;
	unsigned shape_len =  len - ShapeConfigures::instance.get_shape(0).length_ + 1;

	// cout<<"shape len = "<<shape_len<<", current mask = "<<current_mask<<endl;

	while(i < shape_len) 
	{
		if(len-i > 32)
		{
			mask |= reduced_match32(query+32,subject+32,len-i-32) << 32;
		}
			
		for(unsigned j=0;j<32 && i<shape_len;++j) 
		{
			assert(&subject[j] >= ReferenceSeqs<_val>::data_->data(0) && &subject[j] <= ReferenceSeqs<_val>::data_->data(ReferenceSeqs<_val>::data_->raw_len()-1));
			for(unsigned k=0;k<sid;++k)
			{
				if(previous_shape_collision<_val,_pos>(mask, ShapeConfigures::instance.get_shape(k).mask_, &subject[j], k))
				{
					return false;
				}
			}
			// bool b = shape_collision_left<_val,_pos>(mask, current_mask, &subject[j], sid, chunked);

			// cout<<"left"<<", shape_collision_left = "<<b<<endl;
			// if (sequence_type() == amino_acid)
			// {
				if(i < seed_offset && shape_collision_left<_val,_pos>(mask, current_mask, &subject[j], sid, chunked))
				{
						return false;
				}
			// }
			// bool a = shape_collision_right<_val,_pos>(mask, current_mask, &subject[j], sid);
			// cout<<"right"<<", shape_collision_right = "<<a<<endl;
			if(chunked && i > seed_offset && shape_collision_right<_val,_pos>(mask, current_mask, &subject[j], sid))
			{

					return false;
			}
				
			++i;
			mask >>= 1;
		}
		query += 32;
		subject += 32;
	}
	return true;
}

#endif /* COLLISION_H_ */
