
#ifndef STATISTICS_H_
#define STATISTICS_H_

typedef uint64_t stat_type;

class Statistics
{
	public:
	enum value { SEED_HITS, TENTATIVE_MATCHES0, TENTATIVE_MATCHES1, TENTATIVE_MATCHES2, TENTATIVE_MATCHES3, MATCHES, ALIGNED, GAPPED, DUPLICATES,
		GAPPED_HITS, QUERY_SEEDS, QUERY_SEEDS_HIT, REF_SEEDS, REF_SEEDS_HIT, QUERY_SIZE, REF_SIZE, OUT_HITS, OUT_MATCHES, COLLISION_LOOKUPS, QCOV, BIAS_ERRORS, SCORE_TOTAL, ALIGNED_QLEN, COUNT };

	Statistics()
	{ memset(data_, 0, sizeof(data_)); }

	Statistics& operator+=(const Statistics &rhs)
	{
		tthread::lock_guard<tthread::mutex> lock (mtx_);
		for(unsigned i=0;i<COUNT;++i)
			data_[i] += rhs.data_[i];
		return *this;
	}

	void inc(const value v, stat_type n = 1lu)
	{ data_[v] += n; }

	stat_type get(const value v) const
	{ return data_[v]; }

	void print() const
	{
		cout << "Seed hits = " << data_[SEED_HITS] << endl;
		// cout << "Pre_Filter = " << data_[TENTATIVE_MATCHES1] << endl;
		cout << "Matches = " << data_[OUT_MATCHES] << endl;
	}

	stat_type data_[COUNT];
	tthread::mutex mtx_;

} statistics;

#endif /* STATISTICS_H_ */
