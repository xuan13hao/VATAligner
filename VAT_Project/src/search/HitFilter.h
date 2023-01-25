
#ifndef HIT_FILTER_H_
#define HIT_FILTER_H_

#include <vector>
#include <boost/thread/tss.hpp>
#include "trace_pt_buffer.h"
#include "../dp/smith_waterman.h"
#include "../basic/sequence.h"

using std::vector;
using boost::thread_specific_ptr;

template<typename _val, typename _locr, typename _locq, typename _locl>
class HitFilter
{
public:
	HitFilter(Statistics &stats,
			   _locq q_pos,
			   typename Trace_pt_buffer<_locr,_locl>::Iterator &out):
		q_num_ (std::numeric_limits<unsigned>::max()),
		seed_offset_ (std::numeric_limits<unsigned>::max()),
		stats_ (stats),
		q_pos_ (q_pos),
		out_ (out),
		subjects_ (&s2)
		//subjects_ (subjects_ptr)
	{ subjects_->clear(); }

	void push(_locr subject, int score)
	{
		//minimum score to keep a tentative alignment default (0)
		if(score >= VATParameters::min_hit_score)
			push_hit(subject);
		else
			subjects_->push_back(ReferenceSeqs<_val>::data_->fixed_window_infix(subject+VATConsts::seed_anchor));

		
	}

	void finish()
	{
		if(subjects_->size() == 0)
			return;
		unsigned left;
		sequence<const _val> query (QuerySeqs<_val>::data_->window_infix(q_pos_ + VATConsts::seed_anchor, left));
		smith_waterman(query,
				*subjects_,
				VATParameters::hit_band,
				left,
				VATParameters::gap_open + VATParameters::gap_extend,
				VATParameters::gap_extend,
				VATParameters::min_hit_score,
				*this,
				uint8_t(),
				stats_);
	}

	void push_hit(_locr subject)
	{
		if(q_num_ == std::numeric_limits<unsigned>::max()) 
		{
			std::pair<size_t,size_t> l (QuerySeqs<_val>::data_->local_position(q_pos_));
			q_num_ = l.first;
			seed_offset_ = l.second;
		}
		
		assert(subject < ReferenceSeqs<_val>::get().raw_len());
		//seed offset =  suject_end_position - subject_start_point 
		// cout<<"q_num = "<<q_num_<<", subject = "<<subject<<",seed offset = "<<seed_offset_<<endl;
		out_.push(Hits<_locr,_locl> (q_num_, subject, seed_offset_));
		// cout<<"q_num = "<<q_num_<<", subject = "<<subject<<",seed offset = "<<seed_offset_<<endl;
		stats_.inc(Statistics::TENTATIVE_MATCHES3);
	}

	void operator()(int i, const sequence<const _val> &seq, int score)
	{ 
		push_hit(ReferenceSeqs<_val>::data_->position(seq.data()+VATParameters::window-VATConsts::seed_anchor)); 
	
		stats_.inc(Statistics::GAPPED_HITS); 
	}

private:

	unsigned q_num_, seed_offset_;
	Statistics  &stats_;
	_locq q_pos_;
	typename Trace_pt_buffer<_locr,_locl>::Iterator &out_;
	//Tls<vector<sequence<const _val> > > subjects_;
	vector<sequence<const _val> > s2;
	vector<sequence<const _val> >* subjects_;

	//static thread_specific_ptr<vector<sequence<const _val> > > subjects_ptr;

};

/*template<typename _val, typename _locr, typename _locq, typename _locl>
thread_specific_ptr<vector<sequence<const _val> > > hit_filter<_val,_locr,_locq,_locl>::subjects_ptr;*/

#endif /* HIT_FILTER_H_*/
