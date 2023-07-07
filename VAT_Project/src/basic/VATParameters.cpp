

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "VATParameters.h"
#include "AlphabetType.h"

namespace VATParameters {

string		input_ref_file;
uint32_t	threads_;
string		database;
string		query_file;
uint32_t	merge_seq_treshold;
uint32_t	block_size;
uint32_t	hit_cap;
int			min_ungapped_raw_score;
uint32_t	shapes;
uint32_t	index_mode;
uint64_t	max_alignments;
string		output_file;
string		match_file1;
string		match_file2;
int			padding;
uint32_t	output_threads;
uint32_t	compression;
uint32_t	lowmem;
double		chunk_size;
unsigned	min_identities;
unsigned	min_identities2;
int			xdrop;
unsigned	window;
int			min_hit_score;
int			hit_band;
unsigned	min_compressed_identities;
int			min_seed_score;
unsigned	seed_signatures;
double		min_bit_score;
unsigned	run_len;
bool		alignment_traceback;
double		max_seed_freq;
string		tmpdir;
bool		long_mode;
int			gapped_xdrop;
double		max_evalue;
string		sam_output;
string		kegg_file;
int			gap_open;
int			gap_extend;
string		matrix;
string		seg;
bool		verbose;
bool		debug_log;
bool		have_ssse3;
bool		salltitles;
int			reward;
int			penalty;
string		db_type;
double		min_id;
unsigned	compress_temp;
double		toppercent;
string		daa_file;
string		output_format;
bool		forwardonly;
bool		forward_only;
unsigned	fetch_size;
bool		single_domain;
size_t		db_size;
bool		spilce;
bool		chimera;
bool		whole_genome;
bool		circ;
int		seed_len;

Aligner_mode aligner_mode;
Command algn_type;

template<typename _val>
void set_options(double block_size)
{
	if(aligner_mode == accuracy_model) {
		set_option(seed_signatures, 2u);
		set_option(hit_cap, 1024u);
		set_option(index_mode, 3u);
		lowmem = std::max(lowmem, 4u);
	} else
	if(aligner_mode == long_model) {
		set_option(seed_signatures, 1u);
		set_option(index_mode, 2u);
		//lowmem = std::max(lowmem, 4u);
	} else if (aligner_mode == short_model) {
		set_option(seed_signatures, 1u);
		set_option(index_mode, 1u);
	}

	set_option(chunk_size, block_size);

}

string get_temp_file()
{
	if(strlen(getenv("TMPDIR")) > 0)
		return string(getenv("TMPDIR")) + "/VAT.tmp";
	else {
		std::cerr << "Warning: TMPDIR environment variable not set - using output directory for temporary storage.\n";
		return output_file + ".tmp";
	}
}


template<typename _val>
unsigned read_padding(size_t len)
{
	if(padding == 0) {
		if(len<=255)
			return 10;
		else
			return 32;
	} else
		return padding;
}

template<>
unsigned read_padding<Protein>(size_t len)
{
	if(padding == 0) {
		if(len<=35)
			return 5;
		else if(len<=55)
			return 16;
		else
			return 32;
	} else
		return padding;
}
template<>
unsigned read_padding<DNA>(size_t len)
{
	if(padding == 0) {
		if(len<=35)
			return 5;
		else if(len<=55)
			return 16;
		else
			return 32;
	} else
		return padding;
}


template void set_options<Protein>(double block_size);
template void set_options<DNA>(double block_size);
template unsigned read_padding<DNA>(size_t len);
template unsigned read_padding<Protein>(size_t len);

bool mem_buffered()
{ return tmpdir == "/dev/shm"; }

}
