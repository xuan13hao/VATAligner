
#ifndef ALIGN_H_
#define ALIGN_H_

#include "../data/reference.h"
#include "../basic/statistics.h"
#include "../basic/score_matrix.h"
#include "../basic/shape_config.h"
#include "../search/sse_dist.h"
#include "../search/collision.h"
#include "../search/hit_filter.h"
#include "../search/align_ungapped.h"

template<typename _val, typename _locr, typename _locq, typename _locl>
void align(const _locq q_pos,
	  const _val *query,
	  _locr s,
	  Statistics &stats,
	  const unsigned sid,
	  hit_filter<_val,_locr,_locq,_locl> &hf)
{
	stats.inc(Statistics::TENTATIVE_MATCHES0);
	const _val* subject = ref_seqs<_val>::data_->data(s);

	if(fast_match(query, subject) < program_options::min_identities)
		return;

	stats.inc(Statistics::TENTATIVE_MATCHES1);

	unsigned delta, len;
	int score;
	if((score = xdrop_ungapped<_val,_locr,_locq>(query, subject, shape_config::get().get_shape(sid).length_, delta, len)) < program_options::min_ungapped_raw_score)
		return;

	if(!is_primary_hit<_val,_locr>(query-delta, subject-delta, delta, sid, len))
		return;

	stats.inc(Statistics::TENTATIVE_MATCHES2);
	hf.push(s, score);
}

#endif /* ALIGN_H_ */
