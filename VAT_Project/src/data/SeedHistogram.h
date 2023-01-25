

#ifndef SEED_HISTOGRAM_H_
#define SEED_HISTOGRAM_H_

#include "../basic/SeedPartition.h"
#include "SequenceSet.h"
#include "../basic/ShapeParameter.h"
#include "../tools/thread.h"

using std::vector;

void encode_zero_rle(const int32_t *data, size_t len, OutputStreamer &out)
{
	const int32_t *p = data, *end = data + len;
	int32_t n = 0;
	while(p < end) {
		while(p < end && *p == 0) {
			--n;
			++p;
		}
		if(n < 0) {
			out.write(&n, 1);
			n = 0;
		}
		if(p < end) {
			out.write(p, 1);
			++p;
		}
	}
}

void decode_zero_rle(int32_t *data, size_t len, Input_stream &in)
{
	size_t n = 0;
	while(n < len) {
		int32_t x;
		in.read(&x, 1);
		if(x >= 0) {
			*(data++) = x;
			++n;
		} else {
			for(int32_t i=0;i>x;--i) {
				*(data++) = 0;
				++n;
			}
		}
	}
}

typedef int32_t ShapeHistogram[VATConsts::seqp][VATConsts::seedp];

struct seedp_range
{
	seedp_range():
		begin_ (0),
		end_ (0)
	{ }
	seedp_range(unsigned begin, unsigned end):
		begin_ (begin),
		end_ (end)
	{ }
	bool contains(unsigned i) const
	{ return i >= begin_ && i < end_; }
	unsigned begin() const
	{ return begin_; }
	unsigned end() const
	{ return end_; }
	bool lower(unsigned i) const
	{ return i < begin_; }
	bool lower_or_equal(unsigned i) const
	{ return i < end_; }
private:
	unsigned begin_, end_;
} current_range;
//shape histogram partition size 
size_t partition_size(const ShapeHistogram &hst, unsigned p)
{
	size_t s (0);
	for(unsigned i=0;i<VATConsts::seqp;++i)
		s += hst[i][p];
	return s;
}

size_t hst_size(const ShapeHistogram &hst, const seedp_range &range)
{
	size_t s (0);
	for(unsigned i=range.begin();i<range.end();++i)
		s += partition_size(hst, i);
	return s;
}

class SeedHistogram
{
	public:
	SeedHistogram()
	{ }

	template<typename _val>
	SeedHistogram(const SequenceSet<_val> &seqs, const _val&)
	{
		memset(data_, 0, sizeof(data_));
		Build_context<_val> context (seqs, *this);
		launch_scheduled_thread_pool(context, VATConsts::seqp, VATParameters::threads());
	}

	const ShapeHistogram& get(unsigned index_mode, unsigned sid) const
	{ return data_[index_mode-1][sid]; }

	size_t max_chunk_size() const
	{
		size_t max (0);
		::partition p (VATConsts::seedp, VATParameters::lowmem);
		for(unsigned shape=0;shape < ShapeConfigures::get().count();++shape)
			for(unsigned chunk=0;chunk < p.parts; ++chunk)
				max = std::max(max, hst_size(data_[VATParameters::index_mode-1][shape], seedp_range(p.getMin(chunk), p.getMax(chunk))));
		return max;
	}

	void save(OutputStreamer &out) const
	{
		encode_zero_rle(reinterpret_cast<const int32_t*>(data_), sizeof(data_)/sizeof(int32_t), out);
	}

	void load(Input_stream &in)
	{
		decode_zero_rle(reinterpret_cast<int32_t*>(data_), sizeof(data_)/sizeof(int32_t), in);
	}

private:

	template<typename _val>
	struct Build_context
	{
		Build_context(const SequenceSet<_val> &seqs, SeedHistogram &hst):
			seqs (seqs),
			cfgs (shape_configs<_val>()),
			seq_partition (seqs.partition()),
			hst (hst)
		{ }
		void operator()(unsigned thread_id, unsigned seqp) const
		{ hst.build_seq_partition(seqs, seqp, seq_partition[seqp], seq_partition[seqp+1], cfgs); }
		const SequenceSet<_val> &seqs;
		const vector<ShapeConfigures> cfgs;
		const vector<size_t> seq_partition;
		SeedHistogram &hst;
	};

	template<typename _val>
	void build_seq_partition(const SequenceSet<_val> &seqs,
			const unsigned seqp,
			const size_t begin,
			const size_t end,
			const vector<ShapeConfigures> &cfgs)
	{
		// int count = 0;
		assert(seqp < VATConsts::seqp);
		uint64_t key;
		for(size_t i=begin;i<end;++i) {

			assert(i < seqs.get_length());
			const sequence<const _val> seq = seqs[i];
			if(seq.length() < VATConsts::min_shape_len) continue;
			for(unsigned j=0;j<seq.length()+1-VATConsts::min_shape_len; ++j)
				for(vector<ShapeConfigures>::const_iterator cfg = cfgs.begin(); cfg != cfgs.end(); ++cfg) {
					assert(cfg->mode() < VATConsts::index_modes);
					assert(cfg->count() <= VATConsts::max_shapes);
					for(unsigned k=0;k<cfg->count(); ++k)
						if(j+cfg->get_shape(k).length_ < seq.length()+1 && cfg->get_shape(k).set_seed(key, &seq[j]))
						{
							// count++;
							++data_[cfg->mode()][k][seqp][seed_partition(key)];
							// cout<<"index = "<<count<<", data = "<<****data_<<", model = "<<cfg->mode()<<",k = "<<k<<",seqp = "<<seqp<<", key = "<<seed_partition(key)<<endl;
						}

				}

		}
	}

	template<typename _val>
	static vector<ShapeConfigures> shape_configs()
	{
		vector<ShapeConfigures> v;
		if(VATParameters::algn_type == VATParameters::makevatdb) {
			for(unsigned i=1;i<=VATConsts::index_modes;++i)
				v.push_back(ShapeConfigures (i, _val()));
		} else
			v.push_back(ShapeConfigures (VATParameters::index_mode, _val()));
		return v;
	}

	ShapeHistogram data_[VATConsts::index_modes][VATConsts::max_shapes];

};

#endif /* SEED_HISTOGRAM_H_ */
