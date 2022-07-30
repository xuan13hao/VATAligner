#ifndef __ALIGN_FECTCHER_H__
#define __ALIGN_FECTCHER_H__
#include <memory>
#include <vector>
#include <stdint.h>
#include "../basic/sequence.h"
#include "../basic/Hits.h"
#include "../tools/Queue.h"
#include "../tools/text_buffer.h"
#include "../tools/flat_array.h"
using std::vector;
using std::unique_ptr;
using std::list;
using std::array;
/*
struct SeedHit {
	int i, j;
	unsigned frame;
};

struct Parameters
{
	Parameters(uint64_t db_seqs, uint64_t db_letters):
		db_seqs(db_seqs),
		db_letters(db_letters)
	{}
	const uint64_t db_seqs, db_letters;
};


struct Match 
{
	Match(size_t target_block_id, bool outranked):
		target_block_id(target_block_id),
		filter_score(0),
		outranked(outranked)
	{}
	void add_hit(std::list<Hsp> &list, std::list<Hsp>::iterator it) {
		hsp.splice(hsp.end(), list, it);
	}
	bool operator<(const Match &m) const {
		return filter_score > m.filter_score || (filter_score == m.filter_score && target_block_id < m.target_block_id);
	}
	Match(size_t target_block_id, bool outranked, std::array<std::list<Hsp>, MAX_CONTEXT> &hsp);
	void inner_culling(int source_query_len);
	void max_hsp_culling();
	void apply_filters(int source_query_len, const char *query_title);
	size_t target_block_id;
	int filter_score;
	bool outranked;
	std::list<Hsp> hsp;
};

std::vector<Match> extend(const Parameters &params, size_t query_id, Trace_pt_list::iterator begin, Trace_pt_list::iterator end, Statistics &stat, int flags);
TextBuffer* generate_output(vector<Match> &targets, size_t query_block_id, const Parameters &parameters);

}

std::vector<Match> extend(const Parameters &params, size_t query_id, Trace_pt_list::iterator begin, Trace_pt_list::iterator end,  Statistics &stat, int flags);


namespace Extension {



WorkTarget ungapped_stage(const SeedHit *begin, const SeedHit *end, const sequence *query_seq, const Bias_correction *query_cb, size_t block_id) {
	array<vector<Diagonal_segment>, MAX_CONTEXT> diagonal_segments;
	WorkTarget target(block_id, ref_seqs::get()[block_id]);
	for (const SeedHit *hit = begin; hit < end; ++hit) {
		const auto d = xdrop_ungapped(query_seq[hit->frame], target.seq, hit->i, hit->j);
		if(d.score >= config.min_ungapped_raw_score) diagonal_segments[hit->frame].push_back(d);
	}
	for (unsigned frame = 0; frame < align_mode.query_contexts; ++frame) {
		if (diagonal_segments[frame].empty())
			continue;
		std::stable_sort(diagonal_segments[frame].begin(), diagonal_segments[frame].end(), Diagonal_segment::cmp_diag);
		pair<int, list<Hsp_traits>> hsp = greedy_align(query_seq[frame], query_cb[frame], target.seq, diagonal_segments[frame].begin(), diagonal_segments[frame].end(), config.log_extend, frame);
		target.filter_score = std::max(target.filter_score, hsp.first);
		target.hsp[frame] = std::move(hsp.second);
		target.hsp[frame].sort(Hsp_traits::cmp_diag);
	}
	return target;
}

void ungapped_stage_worker(size_t i, size_t thread_id, const sequence *query_seq, const Bias_correction *query_cb, const FlatArray<SeedHit> *seed_hits, size_t *target_block_ids, vector<WorkTarget> *out, mutex *mtx) {
	WorkTarget target = ungapped_stage(seed_hits->begin(i), seed_hits->end(i), query_seq, query_cb, target_block_ids[i]);
	{
		std::lock_guard<mutex> guard(*mtx);
		out->push_back(std::move(target));
	}
}

vector<WorkTarget> ungapped_stage(const sequence *query_seq, const Bias_correction *query_cb, Trace_pt_list::iterator begin, Trace_pt_list::iterator end, int flags) {
	vector<WorkTarget> targets;
	if (begin >= end)
	{
		return targets;
	}
		
	std::sort(begin, end, hit::cmp_subject);
	// cout<<"begin = "<<begin.<<endl;
	size_t target = SIZE_MAX;
	thread_local FlatArray<SeedHit> hits;
	thread_local vector<size_t> target_block_ids;
	hits.clear();
	target_block_ids.clear();
	for (Trace_pt_list::iterator i = begin; i < end; ++i) 
	{
		std::pair<size_t, size_t> l = ref_seqs::data_->local_position(i->subject_);
		if (l.first != target) 
		{
			hits.next();
			target = l.first;
			target_block_ids.push_back(target);
		}
		cout<<" seed_offset_= "<<i->seed_offset_<<", s = "<<l.second<<", q = "<<i->query_<<endl;
		hits.push_back({ (int)i->seed_offset_, (int)l.second, i->query_ % align_mode.query_contexts });
	}

	if (flags & TARGET_PARALLEL) {
		mutex mtx;
		Util::Parallel::scheduled_thread_pool_auto(config.threads_, hits.size(), ungapped_stage_worker, query_seq, query_cb, &hits, target_block_ids.data(), &targets, &mtx);
	}
	else
		for (size_t i = 0; i < hits.size(); ++i)
			targets.push_back(ungapped_stage(hits.begin(i), hits.end(i), query_seq, query_cb, target_block_ids[i]));

	return targets;
}

}
*/

template<typename _locr, typename _locl>
class Align_fetcher
{
	public:
	static void init(size_t qbegin, size_t qend, typename vector<Hits<_locr,_locl> >::iterator begin, typename vector<Hits<_locr,_locl> >::iterator end)
	{
		it_ = begin;
		end_ = end;
		queue_ = unique_ptr<Queue>(new Queue(qbegin, qend));
	}
	bool operator()(size_t query)
	{
		const unsigned q = (unsigned)query,
			c = 1;
		begin = it_;
		while (it_ < end_ && it_->query_ / c == q)
			++it_;
		end = it_;
		this->query = query;
		target_parallel = (end - begin > 0);
		// target_parallel = (end - begin > 0) && (config.frame_shift == 0 || (config.toppercent < 100 && config.query_range_culling));
		return target_parallel;
	}
	bool get()
	{
		return queue_->get(*this) != Queue::end;
	}
	void release() {
		if (target_parallel)
			queue_->release();
	}
	size_t query;
	typename vector<Hits<_locr,_locl> >::iterator begin, end;
	bool target_parallel;
private:	
	static typename vector<Hits<_locr,_locl> >::iterator it_, end_;
	static unique_ptr<Queue> queue_;
};

struct Consumer {
	virtual void consume(const char *ptr, size_t n) = 0;
	virtual void finalize() {}
	virtual ~Consumer() = default;
};

class OutputSink
{
    public:
	OutputSink(size_t begin, Consumer *f) :
		f_(f),
		begin_(begin),
		next_(begin),
		size_(0),
		max_size_(0)
	{}
	void push(size_t n, Text_buffer *buf);
	size_t size() const
	{
		return size_;
	}
	size_t max_size() const
	{
		return max_size_;
	}
	static OutputSink& get()
	{
		return *instance;
	}
	size_t next() const
	{
		return next_;
	}
	size_t begin() const {
		return begin_;
	}
	static std::unique_ptr<OutputSink> instance;
private:
	void flush(Text_buffer *buf);
	std::mutex mtx_;
	Consumer* const f_;
	std::map<size_t, Text_buffer*> backlog_;
	size_t begin_, next_, size_, max_size_;
};
/*
template<typename _val, typename _locr, typename _locl> 
void align_worker(size_t thread_id, const Parameters *params)
{
	Align_fetcher<_locr,_locl> hits;
	while (hits.get()) {
		vector<Extension::Match> matches = Extension::extend(*params, hits.query, hits.begin, hits.end, *metadata, stat, hits.target_parallel ? Extension::TARGET_PARALLEL : 0);
		TextBuffer *buf = Extension::generate_output(matches, hits.query, stat, *metadata, *params);
		if (!matches.empty() && (!config.unaligned.empty() || !config.aligned_file.empty())) {
			query_aligned_mtx.lock();
			query_aligned[hits.query] = true;
			query_aligned_mtx.unlock();
		}
		OutputSink::get().push(hits.query, buf);
		hits.release();
	}
}
*/

template<typename _locr, typename _locl> unique_ptr<Queue> Align_fetcher<_locr,_locl>::queue_;
template<typename _locr, typename _locl> typename vector<Hits<_locr,_locl> >::iterator Align_fetcher<_locr,_locl>::it_;
template<typename _locr, typename _locl> typename vector<Hits<_locr,_locl> >::iterator Align_fetcher<_locr,_locl>::end_;




#endif // __ALIGN_FECTCHER_H__