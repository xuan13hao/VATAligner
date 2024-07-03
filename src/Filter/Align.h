

#ifndef ALIGN_H_
#define ALIGN_H_

#include "../Database/Reference.h"
#include "../Commons/Statistics.h"
#include "../Commons/NuclScoreMatrix.h"
// #include "../basic/ScoreMatrix.h"
#include "../Commons/ShapeParameter.h"
#include "../Filter/sse_dist.h"
#include "../Filter/MainHit.h"
#include "../Filter/HitFilter.h"
#include "../Filter/XdropUngapped.h"

template<typename _val, typename _locr, typename _locq, typename _locl>
void align(const _locq q_pos,
	  const _val *query,
	  _locr s,
	  Statistics &stats,
	  const unsigned sid,
	  HitFilter<_val,_locr,_locq,_locl> &hf)
{
	stats.inc(Statistics::TENTATIVE_MATCHES0);
	//const _val* query = QuerySeqs<_val>::data_->data(q_pos);
	const _val* subject = ReferenceSeqs<_val>::data_->data(s);

	if(fast_match(query, subject) < VATParameters::min_identities)
		return;

	stats.inc(Statistics::TENTATIVE_MATCHES1);

	unsigned delta, len;
	int score;
	if((score = xdrop_ungapped<_val,_locr,_locq>(query, subject, ShapeConfigures::get().get_shape(sid).length_, delta, len)) < VATParameters::min_ungapped_raw_score)
		return;

	// if (sequence_type() == amino_acid)
	// {
	// 	if(!is_primary_hit<_val,_locr>(query-delta, subject-delta, delta, sid, len))
	// 	return;
	// }

	// cout<<"len = "<<len<<", xdrop = "<<score<<", shape id = "<<sid<<", delta = "<<delta<<",shape len = "<<ShapeConfigures::get().get_shape(sid).length_<<endl;

	if(!is_primary_hit<_val,_locr>(query-delta, subject-delta, delta, sid, len))
		return;

	stats.inc(Statistics::TENTATIVE_MATCHES2);
	hf.push(s, score);

}

#endif /* ALIGN_H_ */
