
#ifndef CONST_H_
#define CONST_H_

class VATConsts
{
	public:
	enum {
		build_version = 1,
		build_compatibility = 1,
		daa_version = 0,
		seedp_bits = 10,
		seedp = 1<<seedp_bits,
		max_seed_weight = 32,
		seqp_bits = 8,
		seqp = 1<<seqp_bits,
		max_shapes = 5,
		index_modes = 3,
		min_shape_len = 10,
		max_shape_len = 32,
		seed_anchor = 8
	};

	static const char* version_string;
	static const char* program_name;
	static const char* id_delimiters;

};

const char* VATConsts::version_string = "0.1.01";
const char* VATConsts::program_name = "VAT";
const char* VATConsts::id_delimiters = " \a\b\f\n\r\t\v";

#endif /* CONST_H_ */
