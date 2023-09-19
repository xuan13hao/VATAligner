
#ifndef PACKED_TRANSCRIPT_H_
#define PACKED_TRANSCRIPT_H_

typedef enum { op_match=0, op_insertion=1, op_deletion=2, op_substitution=3 } Edit_operation;

struct Packed_operation
{
	Packed_operation(uint8_t code):
		code (code)
	{ }
	Packed_operation(Edit_operation op, unsigned count):
		code ((op<<6) | count)
	{ }
	template<typename _val>
	Packed_operation(Edit_operation op, _val v):
		code ((op<<6) | (int)v)
	{ }
	operator uint8_t() const
	{ return code; }
	Edit_operation op() const
	{ return (Edit_operation)(code>>6); }
	unsigned count() const
	{ return code&63; }
	template<typename _val>
	_val letter() const
	{ return code&63; }
	static Packed_operation terminator()
	{ return Packed_operation(op_match, 0); }
	uint8_t code;
};

template<typename _val>
struct Combined_operation
{
	Edit_operation op;
	unsigned count;
	_val letter;
};

struct Packed_transcript
{

	template<typename _val>
	struct Const_iterator
	{
		Const_iterator(const Packed_operation *op):
			ptr_ (op)
		{ gather(); }
		bool good() const
		{ return *ptr_ != Packed_operation::terminator(); }
		Const_iterator& operator++()
		{ ++ptr_; gather(); return *this; }
		const Combined_operation<_val>& operator*() const
		{ return op_; }
		const Combined_operation<_val>* operator->() const
		{ return &op_; }
	private:
		void gather()
		{
			if(!good())
				return;
			op_.op = ptr_->op();
			if(op_.op == op_deletion || op_.op == op_substitution) {
				op_.letter = ptr_->letter<_val>();
				op_.count = 1;
			} else {
				op_.count = 0;
				do {
					op_.count += ptr_->count();
					++ptr_;
				} while(good() && ptr_->op() == op_.op);
				--ptr_;
			}
		}
		const Packed_operation *ptr_;
		Combined_operation<_val> op_;
	};

	void read(Buffered_file &f)
	{
		data_.clear();
		uint8_t code;
		do {
			f.read(code);
			data_.push_back(code);
		} while (code != Packed_operation::terminator());
	}

	void read(Binary_buffer::Iterator &it)
	{
		data_.clear();
		uint8_t code;
		do {
			it >> code;
			data_.push_back(code);
		} while (code != Packed_operation::terminator());
	}

	template<typename _val>
	Const_iterator<_val> begin() const
	{ return Const_iterator<_val> (data_.data()); }

	const vector<Packed_operation>& data() const
	{ return data_; }

private:

	vector<Packed_operation> data_;

};

#endif /* PACKED_TRANSCRIPT_H_ */
