/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/


#ifndef SSE_SW_H_
#define SSE_SW_H_

#include "score_profile.h"
#include "dp_matrix.h"

template<typename _score>
inline score_vector<_score> cell_update(const score_vector<_score> &diagonal_cell,
						 const score_vector<_score> &scores,
						 const score_vector<_score> &gap_extension,
						 const score_vector<_score> &gap_open,
						 score_vector<_score> &horizontal_gap,
						 score_vector<_score> &vertical_gap,
						 score_vector<_score> &best,
						 const score_vector<_score> &vbias)
{
	score_vector<_score> current_cell = diagonal_cell + scores;
	current_cell.unbias(vbias);
	current_cell.max(vertical_gap).max(horizontal_gap);
	best.max(current_cell);
	vertical_gap -= gap_extension;
	horizontal_gap -=  gap_extension;
	score_vector<_score> open = current_cell - gap_open;
	vertical_gap.max(open);
	horizontal_gap.max(open);
	return current_cell;
}

template<typename _val, typename _score, typename _callback>
void smith_waterman(const sequence<const _val> &query,
			const vector<sequence<const _val> > &subjects,
			unsigned band,
			unsigned padding,
			int op,
			int ep,
			int filter_score,
			_callback &f,
			const _score&,
			Statistics &stats)
{
	#ifdef SW_ENABLE_DEBUG
	int v[1024][1024];
	#endif

	typedef score_vector<_score> sv;

	unsigned qlen (query.length());
	unsigned slen (subjects[0].length());
	DP_matrix<_score> dp (slen, qlen, band, padding);

	sv open_penalty (static_cast<char>(op));
	sv extend_penalty (static_cast<char>(ep));
	sv vbias (ScoreMatrix::get().bias());
	sequence_stream dseq;
	score_profile<_score> profile;

	typename vector<sequence<const _val> >::const_iterator subject_it (subjects.begin());
	while(subject_it < subjects.end()) 
	{

		const unsigned n_subject (std::min((unsigned)score_traits<_score>::channels, (unsigned)(subjects.end() - subject_it)));
		typename vector<sequence<const _val> >::const_iterator subject_end (subject_it + n_subject);
		sv best;
		dseq.reset();
		dp.clear();

		for(unsigned j=0;j<slen;++j) 
		{
			typename DP_matrix<_score>::Column_iterator it (dp.begin(j));
			sv vgap, hgap, column_best;
			profile.template set_<_val> (dseq.get<_val>(subject_it, subject_end, j, _score()));

			while(!it.at_end()) 
			{
				hgap = it.hgap();
				sv next = cell_update<_score>(it.diag(), profile.get(query[it.row_pos_]), extend_penalty, open_penalty, hgap, vgap, column_best, vbias);
				it.set_hgap(hgap);
				it.set_score(next);
				#ifdef SW_ENABLE_DEBUG
				v[j][it.row_pos_] = next[0];
				#endif
				++it;
			}
			best.max(column_best);
		}

		for(unsigned i=0;i<n_subject;++i)
			if(best[i] >= filter_score)
				f(i + (subject_it - subjects.begin()), *(subject_it + i), best[i]);
		subject_it += score_traits<_score>::channels;
	}

	#ifdef SW_ENABLE_DEBUG
	for(unsigned j=0;j<qlen;++j) {
		for(unsigned i=0;i<subjects[0].length();++i)
			printf("%4i", v[i][j]);
		printf("\n");
	}
	printf("\n");
	#endif
}

#endif /* SSE_SW_H_ */
