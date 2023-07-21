

#ifndef MATCH_H_
#define MATCH_H_
#include <cmath>        // std::abs
#include "sequence.h"
#include "../tools/async_buffer.h"
//include "../tools/util.h"
#include "EditTranscript.h"

enum Strand { FORWARD, REVERSE };

interval normalized_range(unsigned pos, int len, Strand strand)
{
	return strand == FORWARD
			? interval (pos, pos + len)
			: interval (pos + 1 + len, pos + 1);
}

template<typename _locr, typename _locl>
class Hits
{
	public:
	typedef typename packed_sequence_location<_locr>::type packed_loc;

	unsigned	query_;
	packed_loc	subject_;//reference 
	_locl		seed_offset_;//seed offset
	Hits():
		query_ (),
		subject_ (),
		seed_offset_ ()
	{ }
	Hits(unsigned query, _locr subject, _locl seed_offset):
		query_ (query),
		subject_ (subject),
		seed_offset_ (seed_offset)
	{ 
		// cout<<"init hit.."<<endl;
	}
	bool operator<(const Hits &rhs) const
	{ return query_ < rhs.query_; }
	bool blank() const
	{ 
		return subject_ == 0; 
	}
	unsigned operator%(unsigned i) const
	{ 
		return (query_/6) % i; 
	}
	unsigned operator/(unsigned i) const
	{ 
		return (query_/6)/i; 
	}
	int64_t global_diagonal() const
	{ 
		return (int64_t)subject_ - (int64_t)seed_offset_; 
	}
	template<unsigned _d>
	static unsigned query_id(const Hits& x)
	{
 
		return x.query_/_d; 

	}
	template<unsigned _d>
	struct Query_id
	{
		unsigned operator()(const Hits& x) const
		{ return query_id<_d>(x); }
	};
	static bool cmp_subject(const Hits &lhs, const Hits &rhs)
	{ return lhs.subject_ < rhs.subject_; }
	static bool cmp_normalized_subject(const Hits &lhs, const Hits &rhs)
	{
		const uint64_t x = (uint64_t)lhs.subject_ + (uint64_t)rhs.seed_offset_, y = (uint64_t)rhs.subject_ + (uint64_t)lhs.seed_offset_;
		return x < y || (x == y && lhs.seed_offset_ < rhs.seed_offset_);
	}
	friend std::ostream& operator<<(std::ostream &s, const Hits &me)
	{
		s << me.query_ << '\t' << me.subject_ << '\t' << me.seed_offset_ << '\n';
		return s;
	}
} __attribute__((packed));

template<typename _val>
struct local_match
{
	typedef typename vector<local_match>::iterator iterator;
	local_match():
		len_ (),
		query_begin_ (),
		subject_len_ (),
		gap_openings_ (),
		identities_ (),
		mismatches_ (),
		subject_begin_ (),
		score_ (),
		query_len_ (),
		query_anchor_ (0),
		subject_ (0)
	{ }
	local_match(int score):
		len_ (0),
		query_begin_ (0),
		subject_len_ (0),
		gap_openings_ (0),
		identities_ (0),
		mismatches_ (0),
		subject_begin_ (0),
		score_ (score),
		query_len_ (0),
		query_anchor_ (0),
		subject_ (0)
	{ }
	local_match(int query_anchor, const _val *subject):
		len_ (0),
		query_begin_ (0),
		subject_len_ (0),
		gap_openings_ (0),
		identities_ (0),
		mismatches_ (0),
		subject_begin_ (0),
		score_ (0),
		query_len_ (0),
		query_anchor_ (query_anchor),
		subject_ (subject)
	{ }
	local_match(unsigned len, unsigned query_begin, unsigned query_len, unsigned subject_len, unsigned gap_openings, unsigned identities, unsigned mismatches, signed subject_begin, signed score):
		len_ (len),
		query_begin_ (query_begin),
		subject_len_ (subject_len),
		gap_openings_ (gap_openings),
		identities_ (identities),
		mismatches_ (mismatches),
		subject_begin_ (subject_begin),
		score_ (score),
		query_len_ (query_len),
		query_anchor_ (0),
		subject_ (0)
	{ }
	local_match& operator+=(const local_match& rhs)
	{
		add(rhs);
		transcript_right_ = rhs.transcript_right_;
		return *this;
	}
	local_match& operator-=(const local_match& rhs)
	{
		add(rhs);
		query_begin_ = rhs.query_len_;
		subject_begin_ = rhs.subject_len_;
		transcript_left_ = rhs.transcript_right_;
		return *this;
	}
	void add(const local_match &rhs)
	{
		len_ += rhs.len_;
		subject_len_ += rhs.subject_len_;
		gap_openings_ += rhs.gap_openings_;
		identities_ += rhs.identities_;
		mismatches_ += rhs.mismatches_;
		score_ += rhs.score_;
		query_len_ += rhs.query_len_;
	}
	interval query_range(Strand strand) const
	{ return normalized_range(query_begin_, query_len_, strand); }
	interval subject_range() const
	{ return normalized_range(subject_begin_, subject_len_, FORWARD); }
	void print(const sequence<const _val> &query, const sequence<const _val> &subject, const vector<char> &buf) const
	{
		cout << "Score = " << score_ << endl;
		::print(cout, &query[query_begin_], &subject[subject_begin_], transcript_right_, transcript_left_, buf);
	}
	/*friend std::ostream& operator<<(std::ostream &os, const local_match &x)
	{
		os << "(sbj=" << x.subject_range() << " score=" << Scoring<_val>::bitscore(x.score_) << ")";
		return os;
	}*/
	unsigned len_, query_begin_, subject_len_, gap_openings_, identities_, mismatches_;
	signed subject_begin_, score_, query_len_, query_anchor_;
	const _val *subject_;
	Edit_transcript transcript_right_, transcript_left_;
};

template<typename _val>
struct Segment
{
	Segment(int score,
			unsigned frame,
			local_match<_val> *traceback = 0,
			unsigned subject_id = std::numeric_limits<unsigned>::max()):
		score_ (score),
		frame_ (frame),
		traceback_ (traceback),
		subject_id_ (subject_id),
		next_ (0),
		top_score_ (0)
	{ }
	Strand strand() const
	{ return frame_ < 3 ? FORWARD : REVERSE; }
	interval query_range() const
	{ return traceback_->query_range(strand()); }
	interval subject_range() const
	{ return traceback_->subject_range(); }
	bool operator<(const Segment &rhs) const
	{ return top_score_ > rhs.top_score_
			|| (top_score_ == rhs.top_score_
			&& (subject_id_ < rhs.subject_id_ || (subject_id_ == rhs.subject_id_ && (score_ > rhs.score_ || (score_ == rhs.score_ && traceback_->score_ > rhs.traceback_->score_))))); }
	static bool comp_subject(const Segment& lhs, const Segment &rhs)
	{ return lhs.subject_id_ < rhs.subject_id_ || (lhs.subject_id_ == rhs.subject_id_ && lhs.score_ > rhs.score_); }
	struct Subject
	{
		unsigned operator()(const Segment& x) const
		{ return x.subject_id_; }
	};

    // Member function to check strong compatibility with another segment
    bool isStronglyCompatible(const Segment<_val>& other) const
    {
		int MAX_INDEL = 5;
		int MAX_DISTANCE = 50000;
        // Check non-overlapping: Segments should not overlap
        bool nonOverlapping = (query_range().end_ < other.query_range().begin_) ||
                              (other.query_range().end_ < query_range().begin_);

        // Check locality: Segments should be close in query coordinates (e.g., within a certain distance)
        bool locality = static_cast<int>(query_range().end_ - other.query_range().begin_) <= MAX_DISTANCE;

        // Check small INDEL: Compare subject coordinates to check for small INDELs
        bool smallIndel = (static_cast<int>(subject_range().end_ - other.subject_range().begin_) <= MAX_INDEL) ||
                          (static_cast<int>(other.subject_range().end_ - subject_range().begin_) <= MAX_INDEL);

        // Combine the conditions to determine strong compatibility
        return  nonOverlapping && locality && smallIndel;
    }
    bool isChimericMapping(const Segment<_val>& other) const
    {
		int MAX_DISTANCE_IN_REFERENCE = 5;
        // Check if segments correspond to different reference sequences
        bool differentReference = (subject_id_ != other.subject_id_);

        // Check if segments are too far away in a single reference sequence
        bool farAwayInReference = (static_cast<int>(this->subject_range().end_ - other.subject_range().end_) > MAX_DISTANCE_IN_REFERENCE);

        // Determine if segments induce chimeric mapping
        return differentReference || farAwayInReference;
    }
	bool isCompatible(const Segment<_val>& seg1, const Segment<_val>& seg2) {
		// Check if segments are on the same frame (forward or reverse strand)
		if (seg1.frame_ != seg2.frame_)
			return false;

		// Check if the segments are on the same subject sequence (subject_id)
		if (seg1.subject_id_ != seg2.subject_id_)
			return false;

		// Assuming 'local_match' has members query_begin_, query_len_, subject_begin_
		// Check if segments are overlapping and adjacent in the query range
		unsigned overlap = std::min(seg1.traceback_->query_begin_ + seg1.traceback_->query_len_,
									seg2.traceback_->query_begin_ + seg2.traceback_->query_len_)
						- std::max(seg1.traceback_->query_begin_, seg2.traceback_->query_begin_);

		if (overlap > 0) {
			// If overlap is greater than zero, the segments are overlapping
			// Check if the segments are adjacent in the query range
			unsigned gap = seg2.traceback_->query_begin_ - (seg1.traceback_->query_begin_ + seg1.traceback_->query_len_);
			if (gap == 1)
				return true;
		}

		return false;
	}
	int						score_;
	unsigned				frame_;
	local_match<_val>	   *traceback_;
	unsigned				subject_id_;
	Segment				   *next_;
	int						top_score_;
};

#endif /* MATCH_H_ */
