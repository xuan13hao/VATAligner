

#ifndef DAA_WRITE_H_
#define DAA_WRITE_H_

#include "VATFile.h"

class Intermediate_record
{
	public:
	void read(Buffered_file &f)
	{
		f.read(query_id);
		f.read(subject_id);
		f.read(flag);
		f.read_packed(flag & 3, score);
		f.read_packed((flag>>2)&3, query_begin);
		f.read_packed((flag>>4)&3, subject_begin);
		transcript.read(f);
	}
	uint32_t query_id, subject_id, score, query_begin, subject_begin;
	uint8_t flag;
	Packed_transcript transcript;
};

class VATOutput
{
	public:
	VATOutput():
		f_ (VATParameters::daa_file),
		h2_ (ref_header.sequences,
				VATParameters::db_size == 0 ? ref_header.letters : VATParameters::db_size,
				VATParameters::gap_open,
				VATParameters::gap_extend,
				VATParameters::reward,
				VATParameters::penalty,
				ScoreMatrix::get().k(),
				ScoreMatrix::get().lambda(),
				VATParameters::matrix,
				(Align_mode)VATParameters::algn_type)
	{
		VATHeaderOne h1;
		f_.write(&h1, 1);
		h2_.block_type[0] = VATHeaderTwo::alignments;
		h2_.block_type[1] = VATHeaderTwo::ref_names;
		h2_.block_type[2] = VATHeaderTwo::ref_lengths;
		f_.write(&h2_, 1);
	}

	template<typename _val>
	static void write_query_record(Text_buffer &buf, const sequence<const char> &query_name, const sequence<const _val> &query)
	{
		buf.write((uint32_t)0);
		uint32_t l = query.length();
		buf.write(l);
		buf.write_c_str(query_name.c_str(), find_first_of(query_name.c_str(), VATConsts::id_delimiters));
		Packed_sequence s (query);
		uint8_t flags = s.has_n() ? 1 : 0;
		buf.write(flags);
		buf << s.data();
	}

	static void write_record(Text_buffer &buf, const Intermediate_record &r)
	{
		buf.write(r.subject_id).write(r.flag);
		buf.write_packed(r.score);
		buf.write_packed(r.query_begin);
		buf.write_packed(r.subject_begin);
		buf << r.transcript.data();
	}

	template<typename _val>
	static void write_record(Text_buffer &buf,
			const Segment<_val> &match,
			size_t query_source_len,
			const sequence<const _val> &query,
			unsigned query_id,
			const vector<char> &transcript_buf)
	{
		buf.write(ref_map.get<_val>(current_ref_block, match.subject_id_));
		buf.write(get_segment_flag(match));
		buf.write_packed(match.score_);
		buf.write_packed(match.traceback_->query_begin_);
		buf.write_packed(match.traceback_->subject_begin_);
		const unsigned qbegin = query_translated_begin<_val>(match.traceback_->query_begin_, match.frame_, query_source_len, query_translated());

		print_packed(match.traceback_->transcript_right_, match.traceback_->transcript_left_, transcript_buf, buf, query, ReferenceSeqs<_val>::get()[match.subject_id_], qbegin, match.traceback_->subject_begin_);

	}

	void finish()
	{
		uint32_t size = 0;
		f_.write(&size, 1);
		h2_.block_size[0] = f_.tell() - sizeof(VATHeaderOne) - sizeof(VATHeaderTwo);
		h2_.db_seqs_used = ref_map.next_;
		h2_.query_records = statistics.get(Statistics::ALIGNED);

		size_t s = 0;
		for(ptr_vector<string>::const_iterator i = ref_map.name_.begin(); i != ref_map.name_.end(); ++i) {
			f_.write_c_str(*i);
			s += i->length()+1;
		}
		h2_.block_size[1] = s;

		f_.write(ref_map.len_, false);
		h2_.block_size[2] = ref_map.len_.size() * sizeof(uint32_t);

		f_.seek(sizeof(VATHeaderOne));
		f_.write(&h2_, 1);

		f_.close();
	}

	OutputStreamer& stream()
	{ return f_; }

private:

	OutputStreamer f_;
	VATHeaderTwo h2_;

};

#endif /* DAA_WRITE_H_ */
