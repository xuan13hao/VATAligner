

#ifndef OUTPUT_BUFFER_H_
#define OUTPUT_BUFFER_H_

#include "../tools/text_buffer.h"
#include "output_format.h"
#include "VATWrite.h"

unsigned get_length_flag(unsigned x)
{
	if(x <= (unsigned)std::numeric_limits<uint8_t>::max())
		return 0;
	else if(x <= (unsigned)std::numeric_limits<uint16_t>::max())
		return 1;
	else
		return 2;
}

template<typename _val>
unsigned get_rev_flag(unsigned frame)
{ return frame > 0 ? 1 : 0; }

template<>
unsigned get_rev_flag<Protein>(unsigned frame)
{ return frame > 2 ? 1 : 0; }

template<>
unsigned get_rev_flag<DNA>(unsigned frame)
{ return frame > 2 ? 1 : 0; }

template<typename _val>
uint8_t get_segment_flag(const Segment<_val> &match)
{
	unsigned rev = get_rev_flag<_val>(match.frame_);
	return (uint8_t)(get_length_flag(match.score_)
				| (get_length_flag(match.traceback_->query_begin_)<<2)
				| (get_length_flag(match.traceback_->subject_begin_)<<4)
				| rev<<6);
}

template<typename _val>
void write_intermediate_record(Text_buffer &buf,
			const Segment<_val> &match,
			size_t query_source_len,
			const sequence<const _val> &query,
			unsigned query_id,
			const vector<char> &transcript_buf)
{
	buf.write(query_id)
		.write(ref_map.get<_val>(current_ref_block, match.subject_id_))
		.write(get_segment_flag(match))
		.write_packed(match.score_)
		.write_packed(match.traceback_->query_begin_)
		.write_packed(match.traceback_->subject_begin_);
	const unsigned qbegin = query_translated_begin<_val>(match.traceback_->query_begin_, match.frame_, query_source_len, query_translated());
	print_packed(match.traceback_->transcript_right_, match.traceback_->transcript_left_, transcript_buf, buf, query, ReferenceSeqs<_val>::get()[match.subject_id_], qbegin, match.traceback_->subject_begin_);
}

template<typename _val>
class Output_buffer : public Text_buffer
{
	public:
	virtual void print_match(const Segment<_val> &match,
		size_t query_source_len,
		const sequence<const _val> &query,
		unsigned query_id,
		const vector<char> &transcript_buf)
	{ 
		// cout<<"print_match 1"<<endl;
		VATOutput::write_record(*this, match, query_source_len, query, query_id, transcript_buf); 
		// cout<<"print_match 2"<<endl;
	}
	virtual void write_query_record(unsigned query_id)
	{
		query_begin_ = this->size();
		if(query_translated()){
			VATOutput::write_query_record(*this, query_ids::get()[query_id], query_source_seqs::get()[query_id]);
			// cout<<"write_query_record"<<endl;
			}
		else
			VATOutput::write_query_record(*this, query_ids::get()[query_id], QuerySeqs<_val>::get()[query_id]);
	}
	virtual void finish_query_record()
	{ *(uint32_t*)(this->data_+query_begin_) = this->size() - query_begin_ - sizeof(uint32_t); }
	virtual ~Output_buffer()
	{ }
private:
	size_t query_begin_;
};

template<typename _val>
class Temp_output_buffer : public Output_buffer<_val>
{
	public:
	virtual void print_match(const Segment<_val> &match,
				size_t query_source_len,
				const sequence<const _val> &query,
				unsigned query_id,
				const vector<char> &transcript_buf)
	{ write_intermediate_record(*this, match, query_source_len, query, query_id, transcript_buf); }
	virtual void write_query_record(unsigned query_id)
	{ }
	virtual void finish_query_record()
	{ }
	virtual ~Temp_output_buffer()
	{ }
};

#endif /* OUTPUT_BUFFER_H_ */
