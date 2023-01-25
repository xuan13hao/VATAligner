

#ifndef COMPLEXITY_FILTER_H_
#define COMPLEXITY_FILTER_H_

#include "../algo/blast/core/blast_seg.h"
#include "../algo/blast/core/blast_filter.h"
#include "../basic/value.h"

template<class _val>
struct Complexity_filter
{
	unsigned filter(vector<_val> &seq) const
	{ return 0; }
	static const Complexity_filter& get()
	{ return instance; }
	void run(AlphabetSet<_val> &seqs) const
	{ }
	static const Complexity_filter instance;
};

template<>
struct Complexity_filter<Protein>
{

	Complexity_filter()
	{ blast_seg_ = SegParametersNewAa(); }

	~Complexity_filter()
	{ SegParametersFree(blast_seg_); }

	unsigned filter(sequence<Protein> seq) const
	{
		BlastSeqLoc *seg_locs;
		SeqBufferSeg ((uint8_t*) seq.data(), seq.length(), 0, blast_seg_, &seg_locs);
		unsigned nMasked = 0;

		if(seg_locs) {
			BlastSeqLoc *l = seg_locs;
			do {
				for(signed i=l->ssr->left;i<=l->ssr->right;i++) {
					nMasked++;
					seq[i] = AlphabetAttributes<Protein>::MASK_CHAR;
				}
			} while((l=l->next) != 0);
			BlastSeqLocFree(seg_locs);
		}
		return nMasked;
	}

	static const Complexity_filter& get()
	{ return instance; }

	void run(AlphabetSet<Protein> &seqs) const
	{
		Filter_context context (seqs, *this);
		launch_scheduled_thread_pool(context, seqs.get_length(), VATParameters::threads());
	}

private:

	struct Filter_context
	{
		Filter_context(AlphabetSet<Protein> &seqs, const Complexity_filter &filter):
			seqs (seqs),
			filter (filter)
		{ }
		void operator()(unsigned thread_id, unsigned i)
		{
			filter.filter(seqs[i]);
		}
		AlphabetSet<Protein> &seqs;
		const Complexity_filter &filter;
	};

	SegParameters *blast_seg_;

	static const Complexity_filter instance;

};

const Complexity_filter<Protein> Complexity_filter<Protein>::instance;
#ifdef NDEBUG
template<> const Complexity_filter<DNA> Complexity_filter<DNA>::instance;
#else
template<typename _val> const Complexity_filter<_val> Complexity_filter<_val>::instance;
#endif


#endif /* COMPLEXITY_FILTER_H_ */
