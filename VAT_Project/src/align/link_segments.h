
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
	// cout<<"-------------------"<<endl;
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

template<typename _val>
void link_Chimeric_segments(vector<Segment<_val> > &hsp_list)
{
    hsp_list = hsp_list;
}


template<typename _val>
void link_wgs_segments(vector<Segment<_val>>& hsp_list)
{
    // Define a vector to store the detected synteny blocks
    vector<vector<Segment<_val>>> syntenyBlocks;

    // Sort the hsp_list based on the comp_subject comparison function
    std::sort(hsp_list.begin(), hsp_list.end(), Segment<_val>::comp_subject);

    // Iterate over the hsp_list to detect synteny blocks
    for (const Segment<_val>& segment : hsp_list)
    {
        bool addedToBlock = false;

        // Try to add the segment to an existing synteny block
        for (vector<Segment<_val>>& block : syntenyBlocks)
        {
            if (segment.isStronglyCompatible(block.front()))
            {
                block.push_back(segment);
                addedToBlock = true;
                break;
            }
        }

        // If the segment didn't fit into any existing block, create a new block
        if (!addedToBlock)
        {
            syntenyBlocks.push_back({ segment });
        }
    }

    // Order the synteny blocks based on a specific criteria (e.g., top_score_)
    std::sort(syntenyBlocks.begin(), syntenyBlocks.end(),
              [](const vector<Segment<_val>>& block1, const vector<Segment<_val>>& block2) {
                  return block1.front().top_score_ > block2.front().top_score_;
              });

    // Clear the hsp_list vector
    hsp_list.clear();

    // Concatenate the ordered synteny blocks into the hsp_list vector
    for (const vector<Segment<_val>>& block : syntenyBlocks)
    {
        hsp_list.insert(hsp_list.end(), block.begin(), block.end());
    }
}

template<typename _val>
void findchimericSegments(std::vector<Segment<_val>>& segments, int maxDistance) {
    std::vector<int> dp(segments.size());
    std::vector<int> prev(segments.size(), -1);
    std::vector<int> maxLen(segments.size());
    std::vector<int> maxIdx(segments.size());

    int bestScore = 0;
    int bestIdx = -1;

    for (int i = 0; i < segments.size(); ++i) {
        dp[i] = segments[i].score_;
        maxLen[i] = segments[i].score_;
        maxIdx[i] = i;

        for (int j = 0; j < i; ++j) {
            if (!segments[i].isChimericMapping(segments[j]))
                continue;

            int distance = std::abs(static_cast<int>(segments[i].subject_range().begin_) - static_cast<int>(segments[j].subject_range().end_));
            if (distance > maxDistance)
                continue;

            int score = dp[j] + segments[i].score_;
            if (score > dp[i]) {
                dp[i] = score;
                prev[i] = j;
            }
        }

        if (dp[i] > bestScore) {
            bestScore = dp[i];
            bestIdx = i;
        }
    }

    std::vector<std::vector<Segment<_val>>> chainedSegments;
    while (bestIdx >= 0) {
        chainedSegments.push_back({ segments[bestIdx] });
        bestIdx = prev[bestIdx];
    }

    std::reverse(chainedSegments.begin(), chainedSegments.end());
}



#endif /* LINK_SEGMENTS_H_ */
