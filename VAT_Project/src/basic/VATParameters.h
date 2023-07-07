

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <stdint.h>
#include <string>

using std::string;

namespace VATParameters
{

	extern string	input_ref_file;
	extern uint32_t	threads_;
	extern string	database;
	extern string	query_file;
	extern uint32_t	merge_seq_treshold;
	extern uint32_t	block_size;
	extern uint32_t	hit_cap;
	extern int		min_ungapped_raw_score;
	extern uint32_t shapes;
	extern uint32_t	index_mode;
	extern uint64_t	max_alignments;
	extern string	match_file1;
	extern string	match_file2;
	extern int		padding;
	extern uint32_t	output_threads;
	extern uint32_t compression;
	extern uint32_t lowmem;
	extern double	chunk_size;
	extern unsigned min_identities;
	extern unsigned min_identities2;
	extern int		xdrop;
	extern unsigned window;
	extern int		min_hit_score;
	extern int		hit_band;
	extern unsigned	min_compressed_identities;
	extern int		min_seed_score;
	extern unsigned	seed_signatures;
	extern double	min_bit_score;
	extern unsigned	run_len;
	extern bool		alignment_traceback;
	extern double	max_seed_freq;
	extern string	tmpdir;
	extern bool		long_mode;
	extern int		gapped_xdrop;
	extern double	max_evalue;
	extern string	kegg_file;
	extern int		gap_open;
	extern int		gap_extend;
	extern string	matrix;
	extern string	seg;
	extern bool		verbose;
	extern bool		debug_log;
	extern bool		have_ssse3;
	extern bool		salltitles;
	extern int		reward;
	extern int		penalty;
	extern string	db_type;
	extern double	min_id;
	extern unsigned	compress_temp;
	extern double	toppercent;
	extern string	daa_file;
	extern string	output_format;
	extern string	output_file;
	extern bool		forwardonly;
	extern bool		forward_only;
	extern unsigned fetch_size;
	extern bool		single_domain;
	extern size_t	db_size;
	extern bool		spilce;
	extern bool		chimera;
	extern bool		whole_genome;
	extern bool		circ;
	extern int		seed_len;

	typedef enum { short_model=0, long_model=1,accuracy_model=3 } Aligner_mode;
	extern Aligner_mode aligner_mode;
	typedef enum { invalid=0, makevatdb=1, protein=2, blastx=3, dna=4, view=5 } Command;
	extern Command algn_type;

	inline uint32_t threads()
	{
		return std::max(threads_, 1U);
	}

	template<typename _t>
	inline void set_option(_t& option, _t value)
	{
		if(option == 0)
			option = value;
	}

	template<typename _val>
	void set_options(double block_size);
	template<typename _val>
	unsigned read_padding(size_t len);
	string get_temp_file();
	bool mem_buffered();

	inline unsigned get_run_len(unsigned length)
	{
		if(run_len == 0) {
			if(length < 30)
				return 1;
			else if(length < 100)
				return 20;
			else
				return 40;
		} else
			return run_len;
	}

	inline bool output_range(unsigned n_target_seq, int score, int top_score)
	{
		if(toppercent < 100)
			return (1.0-(double)score/top_score)*100 <= toppercent;
		else
			return n_target_seq < max_alignments;
	}

}

#endif /* OPTIONS_H_ */
