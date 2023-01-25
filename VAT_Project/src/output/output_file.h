

#ifndef OUTPUT_FILE_H_
#define OUTPUT_FILE_H_

#include <string>
#include "output_buffer.h"

using std::string;

class Block_output : public Buffered_file
{
public:
	struct Iterator
	{
		unsigned block_;
		bool same_subject_;
		Intermediate_record info_;
		bool operator<(const Iterator &rhs) const
		{ return info_.query_id > rhs.info_.query_id ||
				(info_.query_id == rhs.info_.query_id && (rhs.same_subject_ ||
						(!rhs.same_subject_ && info_.score < rhs.info_.score))); }
	};

	bool next(Iterator &it, unsigned subject, unsigned query)
	{
		if(this->eof())
			return false;
		it.info_.read(*this);
		it.block_ = block_;
		it.same_subject_ = it.info_.subject_id == subject && it.info_.query_id == query;
		return true;
	}

	Block_output(unsigned ref_block, const TempFile &tmp_file):
		Buffered_file (tmp_file),
		block_ (ref_block)
	{ }

private:

	const unsigned block_;

};

#endif /* OUTPUT_FILE_H_ */
