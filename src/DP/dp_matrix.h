

#ifndef DP_MATRIX_H_
#define DP_MATRIX_H_

#include <vector>
#include <boost/thread/tss.hpp>
#include "score_vector.h"

using std::vector;
using boost::thread_specific_ptr;

template<typename _score>
void array_clear(score_vector<_score> *v, unsigned n)
{
	score_vector<_score> *end (v+n);
	while(v < end)
		*(v++) = score_vector<_score> ();
}

template<typename _score>
class DP_matrix
{
	public:
	typedef score_vector<_score> sv;

	struct Column_iterator
	{
		Column_iterator(unsigned column, sv* hgap_front, sv* score_front, unsigned row_pos, unsigned row_end, unsigned delta):
			row_pos_ (row_pos),
			row_end_ (row_end),
			delta_ (delta),
			hgap_ptr_ (hgap_front),
			score_ptr_ (score_front),
			d_ (delta > 0 ? *score_front : sv ())
		{ }
		inline bool at_end() const
		{ return row_pos_ >= row_end_; }
		inline void operator++()
		{ ++row_pos_; ++hgap_ptr_; ++score_ptr_; }
		inline sv hgap() const
		{ return *(hgap_ptr_+delta_); }
		inline sv diag() const
		{ return d_; }
		inline void set_hgap(const sv& x)
		{ *hgap_ptr_ = x; }
		inline void set_score(const sv& x)
		{ d_ = *(score_ptr_+delta_); *score_ptr_ = x; }
		unsigned row_pos_, row_end_, delta_;
		sv *hgap_ptr_, *score_ptr_, d_;
	};

	DP_matrix(unsigned columns, unsigned rows, unsigned band, unsigned padding):
		rows_ (rows),
		band_ (band),
		padding_ (padding),
		scores_ (&s2),
		hgap_ (&h2)
	{
		scores_->resize(2*band+1);
		hgap_->resize(2*band+2);
		hgap_front_ = &hgap_->front();
		score_front_ = &scores_->front();
	}

	inline void clear()
	{
		array_clear(hgap_front_, 2*band_+2);
		array_clear(score_front_, 2*band_+1);
	}

	inline void band_range(unsigned column, unsigned& begin, unsigned& end)
	{
		if(column >= rows_ + padding_) {
			begin = 0;
			end = band_;
		} else if(column >= padding_) {
			unsigned pj (column - padding_);
			unsigned top_delta (pj >= band_ ? 0 : band_ - pj);
			unsigned query_start (pj >= band_ ? pj - band_ : 0);
			unsigned query_end (std::min(pj+band_+1, rows_));
			begin = top_delta;
			end = begin + query_end - query_start;
		} else {
			begin = band_ + 1;
			end = begin + band_;
		}
	}

	inline Column_iterator begin(unsigned column)
	{
		if(column >= rows_ + padding_) {
			return Column_iterator (column, hgap_front_, score_front_, rows_-band_, rows_, 0);
		} else if(column >= padding_) {
			unsigned pj (column - padding_);
			unsigned top_delta (pj >= band_ ? 0 : band_ - pj);
			unsigned query_start (pj >= band_ ? pj - band_ : 0);
			unsigned query_end (std::min(pj+band_+1, rows_));
			return Column_iterator (column, hgap_front_+top_delta, score_front_+top_delta, query_start, query_end, 1);
		} else {
			return Column_iterator (column, hgap_front_+band_+1, score_front_+band_+1, 0, band_, 0);
		}
	}

	inline void sub_all(sv *ptr, const sv *end, const sv& x)
	{
		while(ptr < end)
			*(ptr++) -= x;
	}

	inline sv get_min(const sv *ptr, const sv *end) const
	{
		sv x (*(ptr++));
		while(ptr < end)
			x = x.min(*(ptr++));
		return x;
	}

private:

	static thread_specific_ptr<vector<sv> > scores_ptr;
	static thread_specific_ptr<vector<sv> > hgap_ptr;

	const unsigned rows_, band_, padding_;
	sv *hgap_front_, *score_front_;
	//Tls<vector<sv> > scores_, hgap_;
	vector<sv>* scores_, *hgap_;
	vector<sv> s2,h2;

};

template<typename _score> thread_specific_ptr<vector<score_vector<_score> > > DP_matrix<_score>::scores_ptr;
template<typename _score> thread_specific_ptr<vector<score_vector<_score> > > DP_matrix<_score>::hgap_ptr;

#endif /* DP_MATRIX_H_ */
