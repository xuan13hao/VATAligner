

#ifndef DAA_FILE_H_
#define DAA_FILE_H_

using std::string;

class VATHeaderOne
{
	public:
	VATHeaderOne():
		magic_number (0x3c0e53476d3ee36bllu),
		version (0)
	{ }
	uint64_t magic_number, version;
};

typedef enum { blastp=2, blastx=3, dna=4 } Align_mode;

class VATHeaderTwo
{
	public:
	VATHeaderTwo()
	{ }
	VATHeaderTwo(uint64_t db_seqs,
			uint64_t db_letters,
			int32_t gap_open,
			int32_t gap_extend,
			int32_t reward,
			int32_t penalty,
			double k,
			double lambda,
			const string &score_matrix,
			Align_mode mode):
		vat_build_ (VATConsts::build_version),
		db_seqs (db_seqs),
		db_seqs_used (0),
		db_letters (db_letters),
		flags (0),
		query_records (0),
		mode (mode),
		gap_open (gap_open),
		gap_extend (gap_extend),
		reward (reward),
		penalty (penalty),
		reserved1 (0),
		reserved2 (0),
		reserved3 (0),
		k (k),
		lambda (lambda),
		reserved4 (0),
		reserved5 (0)
	{
		memset(block_type, 0, sizeof(block_type));
		memset(block_size, 0, sizeof(block_size));
		strcpy(this->score_matrix, score_matrix.c_str());
	}
	typedef enum { empty = 0, alignments = 1, ref_names = 2, ref_lengths = 3 } Block_type;
	uint64_t vat_build_, db_seqs, db_seqs_used, db_letters, flags, query_records;
	int32_t mode, gap_open, gap_extend, reward, penalty, reserved1, reserved2, reserved3;
	double k, lambda, reserved4, reserved5;
	char score_matrix[16];
	uint64_t block_size[256];
	char block_type[256];
};

class VATFile
{
public:
	VATFile(const string& file_name):
		f_ (file_name)
	{
		f_.read(&h1_, 1);
		if(h1_.magic_number != VATHeaderOne().magic_number)
			throw std::runtime_error("Input file is not a VAT file.");
		if(h1_.version > VATConsts::daa_version)
			throw std::runtime_error("VAT version requires later version of VAT.");
		f_.read(&h2_, 1);

		if(h2_.block_size[0] == 0)
			throw std::runtime_error("Invalid VAT file. VAT run probably has not completed successfully.");

		f_.seek(sizeof(VATHeaderOne) + sizeof(VATHeaderTwo) + h2_.block_size[0]);
		string s;
		ref_name_.reserve(h2_.db_seqs_used);
		for(uint64_t i=0;i<h2_.db_seqs_used;++i) {
			f_.read_c_str(s);
			ref_name_.push_back(new string(s));
		}
		ref_len_.resize(h2_.db_seqs_used);
		f_.read(ref_len_.data(), h2_.db_seqs_used);

		f_.seek(sizeof(VATHeaderOne) + sizeof(VATHeaderTwo));

	}

	uint64_t diamond_build() const
	{ return h2_.vat_build_; }

	uint64_t db_seqs() const
	{ return h2_.db_seqs; }

	uint64_t db_seqs_used() const
	{ return h2_.db_seqs_used; }

	uint64_t db_letters() const
	{ return h2_.db_letters; }

	const char* score_matrix() const
	{ return h2_.score_matrix; }

	int32_t gap_open_penalty() const
	{ return h2_.gap_open; }

	int32_t gap_extension_penalty() const
	{ return h2_.gap_extend; }

	int32_t match_reward() const
	{ return h2_.reward; }

	int32_t mismatch_penalty() const
	{ return h2_.penalty; }

	uint64_t query_records() const
	{ return h2_.query_records; }

	Align_mode mode() const
	{ return (Align_mode)h2_.mode; }

	const string& ref_name(size_t i) const
	{ return ref_name_[i]; }

	const uint32_t ref_len(size_t i) const
	{ return ref_len_[i]; }

	bool read_query_buffer(Binary_buffer &buf)
	{
		uint32_t size;
		f_.read(&size, 1);
		if(size == 0)
			return false;
		buf.resize(size);
		f_.read(buf.data(), size);
		return true;
	}

private:

	Input_stream f_;
	VATHeaderOne h1_;
	VATHeaderTwo h2_;
	ptr_vector<string> ref_name_;
	vector<uint32_t> ref_len_;

};

#endif /* DAA_FILE_H_ */
