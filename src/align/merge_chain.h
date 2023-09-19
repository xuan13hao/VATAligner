#ifndef __MERGE_CHAIN_H__
#define __MERGE_CHAIN_H__

template<typename _val>
void mergeSegments(Segment<_val>& chain, const Segment<_val>& segment) {
    // Merge the segment into the chain
    chain.score_ += segment.score_;
    
    // Extend the chain's query and subject ranges based on the new segment
    if (segment.traceback_->query_begin_ < chain.traceback_->query_begin_) {
        chain.traceback_->query_len_ += (chain.traceback_->query_begin_ - segment.traceback_->query_begin_);
        chain.traceback_->query_begin_ = segment.traceback_->query_begin_;
    }
    
    unsigned chainQueryEnd = chain.traceback_->query_begin_ + chain.traceback_->query_len_;
    unsigned segmentQueryEnd = segment.traceback_->query_begin_ + segment.traceback_->query_len_;
    if (segmentQueryEnd > chainQueryEnd) {
        chain.traceback_->query_len_ += (segmentQueryEnd - chainQueryEnd);
    }

    if (segment.traceback_->subject_begin_ < chain.traceback_->subject_begin_) {
        chain.traceback_->subject_begin_ = segment.traceback_->subject_begin_;
    }
    
    unsigned chainSubjectEnd = chain.traceback_->subject_begin_ + chain.traceback_->subject_len_; // Assuming local_match has a subject_len_ member
    unsigned segmentSubjectEnd = segment.traceback_->subject_begin_ + segment.traceback_->subject_len_;
    if (segmentSubjectEnd > chainSubjectEnd) {
        chain.traceback_->subject_len_ += (segmentSubjectEnd - chainSubjectEnd);
    }
}

// Function to combine overlapping and adjacent alignment regions into longer chains
template<typename _val>
void combineSegments(std::vector<Segment<_val>>& segments) {
    // Sort the input segments based on their starting positions (query_begin_)
    std::sort(segments.begin(), segments.end(), [](const Segment<_val>& a, const Segment<_val>& b) {
        return a.traceback_->query_begin_ < b.traceback_->query_begin_;
    });

    for (size_t i = 0; i < segments.size(); ++i) {
        for (size_t j = i + 1; j < segments.size(); ++j) {
            if (isCompatible(segments[i], segments[j])) {
                // Merge the segment into the chain
                mergeSegments(segments[i], segments[j]);
                // Remove the merged segment from the vector
                segments.erase(segments.begin() + j);
                --j;
            }
        }
    }
};
#endif // __MERGE_CHAIN_H__