
#ifndef LINK_SEGMENTS_H_
#define LINK_SEGMENTS_H_

#include <vector>
#include "../tools/map.h"

using std::vector;

static const size_t MAX_LINKING_OVERLAP = 10;

template<typename _val>
double link_segments(Segment<_val> &h1, Segment<_val> &h2)
{
	if(h1.strand() == h2.strand() && h2.top_evalue_ == -1) {
		Segment<_val> *p = &h1;
		while(p != 0) {
			if(p->query_range().overlap(h2.query_range())/3 >= 1
				|| p->subject_range().overlap(h2.subject_range()) > MAX_LINKING_OVERLAP)
				return std::numeric_limits<double>::max();
			p = p->next_;
		}
		double ev = h1.evalue_ * h2.evalue_;
		h2.next_ = h1.next_;
		h1.next_ = &h2;
		p = &h1;
		while(p != 0) {
			p->evalue_ = ev;
			p->top_evalue_ = ev;
			p = p->next_;
		}
		return ev;
	}
	return std::numeric_limits<double>::max();
}

template<typename _val>
void link_segments(const typename vector<Segment<_val> >::iterator &begin, const typename vector<Segment<_val> >::iterator &end)
{
	int max_score = begin->score_;
	/*for(typename vector<match<_val> >::iterator i=begin; i<end; ++i)
		if(i->top_evalue_ == -1)
			for(typename vector<match<_val> >::iterator j=i+1; j<end; ++j)
				min_ev = std::min(min_ev, link_segments(*i, *j));*/
	for(typename vector<Segment<_val> >::iterator i=begin; i<end; ++i)
		i->top_score_ = max_score;
}

template<typename _val>
void link_segments(vector<Segment<_val> > &hsp_list)
{
	typedef Map<typename vector<Segment<_val> >::iterator,typename Segment<_val>::Subject> Hsp_map;
	std::sort(hsp_list.begin(), hsp_list.end(), Segment<_val>::comp_subject);
	Hsp_map hsp_map (hsp_list.begin(), hsp_list.end());
	typename Hsp_map::Iterator it = hsp_map.begin();
	while(it.valid()) {
		link_segments<_val>(it.begin(), it.end());
		++it;
	}
}

#endif /* LINK_SEGMENTS_H_ */
