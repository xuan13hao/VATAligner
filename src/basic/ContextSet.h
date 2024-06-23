
#ifndef SETUP_H_
#define SETUP_H_

#include <boost/iostreams/tee.hpp>
#include "VATParameters.h"
#include "../tools/system.h"

using std::cout;
using std::endl;
using std::auto_ptr;
void setup(const string &command, int ac, const char **av)
{
	namespace io = boost::iostreams;
	namespace po = VATParameters;

	auto_append_extension(po::database, ".vatf");
	auto_append_extension(po::daa_file, ".vatr");

	if(po::debug_log) {
		io::tee_filter<io::file_sink> t (io::file_sink("vat.log", std::ios_base::out | std::ios_base::app));
		verbose_stream.push(t);
		log_stream.push(t);
		log_stream.push(cout);
	} else
		log_stream.push(io::null_sink ());
	if(po::verbose || po::debug_log)
		verbose_stream.push(cout);
	else
		verbose_stream.push(io::null_sink ());

	log_stream << "Command line: ";
	for(int i=0;i<ac;++i)
		log_stream << av[i] << ' ';
	log_stream << endl;

	verbose_stream << VATConsts::program_name << " v" << VATConsts::version_string << "." << VATConsts::build_version << endl;
#ifndef NDEBUG
	verbose_stream << "Assertions enabled." << endl;
#endif
	po::set_option(po::threads_, tthread::thread::hardware_concurrency());
	verbose_stream << "#Threads = " << po::threads() << endl;

	if(command == "makevatdb")
		po::algn_type = po::makevatdb;
	else if(command == "blastx")
		po::algn_type = po::blastx;
	else if(command == "protein")
		po::algn_type = po::protein;
	else if(command == "dna")
		po::algn_type = po::dna;
	else if(command == "view")
		po::algn_type = po::view;
	else
		po::algn_type = po::invalid;

	if(sequence_type() == amino_acid) {
		if(po::gap_open == -1)
			po::gap_open = 11;
		if(po::gap_extend == -1)
			po::gap_extend = 1;
		ScoreMatrix::instance = auto_ptr<ScoreMatrix> (new ScoreMatrix(po::matrix,
				po::gap_open,
				po::gap_extend,
				po::reward,
				po::penalty,
				Protein ()));
		// ScoreMatrix::get().print<Protein>();
	} 
	else {
		if(po::gap_open == -1)
			po::gap_open = 5;
		if(po::gap_extend == -1)
			po::gap_extend = 3;

if (po::seed_len == 8)
{
    uint32_t idx = 8 - 8 + 1;
    po::set_option(po::index_mode, idx);
	// cout<<"default = "<<idx<<endl;

}
else if (po::seed_len == 9)
{
    uint32_t idx = 9 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 10)
{
    uint32_t idx = 10 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 11)
{
    uint32_t idx = 11 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 12)
{
    uint32_t idx = 12 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 13)
{
    uint32_t idx = 13 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 14)
{
    uint32_t idx = 14 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 16)
{
    uint32_t idx = 16 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 17)
{
    uint32_t idx = 17 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 18)
{
    uint32_t idx = 18 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 19)
{
    uint32_t idx = 19 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 20)
{
    uint32_t idx = 20 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 21)
{
    uint32_t idx = 21 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 22)
{
    uint32_t idx = 22 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 23)
{
    uint32_t idx = 23 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 24)
{
    uint32_t idx = 24 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 25)
{
    uint32_t idx = 25 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 26)
{
    uint32_t idx = 26 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 27)
{
    uint32_t idx = 27 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 28)
{
    uint32_t idx = 28 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 29)
{
    uint32_t idx = 29 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 30)
{
    uint32_t idx = 30 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else if (po::seed_len == 15)
{
    uint32_t idx = 15 - 8 + 1;
    po::set_option(po::index_mode, idx);
}
else
{
    uint32_t idx = 15 - 8 + 1;
    po::set_option(po::index_mode, idx);
	// cout<<"default = "<<idx<<endl;
}
		if (po::chimera)
		{
			// cout<<"Init chimeric alignment parameters"<<endl;
			// po::set_option(po::hit_cap, 15u);
			// po::set_option(po::min_identities, 14u);
			// po::set_option(po::padding, 8);
			// po::set_option(po::match, 1);
			// po::set_option(po::mismatch, -5);
			// po::set_option(po::lowmem, 1u);
			// cout<<VATParameters::hit_cap<<"\t"<<VATParameters::lowmem<<"\t"<<VATParameters::min_identities<<"\t"<<VATParameters::match<<"\t"<<VATParameters::mismatch<<endl;
			po::lowmem = 1u;
			po::penalty = -4;
			po::gap_extend = 0;
			po::gap_open = 0;
			po::match = 1;
			po::mismatch = -5;
			po::penalty = -5;
			po::reward = 1;
			// po::penalty = -5;
			po::padding = 8;
			po::hit_cap = 15u;
			po::min_identities = 14u;
		}
		else if (po::spilce || po::circ)
		{
			// cout<<"Init spliced alignment parameters"<<endl;
			// po::set_option(po::hit_cap, 15u);
			// po::set_option(po::min_identities, 29u);
			// po::set_option(po::padding, 8);
			// po::set_option(po::match, 5);
			// po::set_option(po::mismatch, -2);
			// po::set_option(po::lowmem, 1u);
			// VATParameters::lowmem = 1u;
			// // VATParameters::penalty = -4;
			// // VATParameters::gap_extend = 0;
			po::gap_extend = 0;
			po::gap_open = 0;
			po::match = 5;
			po::mismatch = -2;
			po::reward = 5;
			po::penalty = -2;
			po::padding = 8;
			po::hit_cap = 15;
			po::min_identities = 28;
		}
		else if (po::whole_genome_sequencing)
		{
			// cout<<"Init whole genome sequencing alignment parameters"<<endl;
			// po::set_option(po::hit_cap, 2u);
			// po::set_option(po::min_identities, 28u);
			// po::set_option(po::padding, 8);
			// // po::set_option(po::match, 5);
			// // po::set_option(po::mismatch, -2);
			// po::set_option(po::lowmem, 1u);
			// VATParameters::lowmem = 4;
			// // VATParameters::penalty = -4;
			// // VATParameters::gap_extend = 0;
			po::lowmem = 1u;
			po::match = 5;
			po::mismatch = -2;
			po::penalty = -2;
			po::reward = 5;
			po::penalty = -2;
			po::padding = 8;
			po::hit_cap = 2;
			po::min_identities = 28;
		}
		else if (po::whole_genome)
		{
			// cout<<"Init whole genome alignment parameters"<<endl;
			// po::set_option(po::hit_cap, 2u);
			// po::set_option(po::min_identities, 20u);
			// po::set_option(po::padding, 8);
			// // po::set_option(po::match, 5);
			// // po::set_option(po::mismatch, -2);
			// po::set_option(po::lowmem, 4u);
			// VATParameters::gap_open = 0;
			// VATParameters::lowmem = 4;
			// VATParameters::penalty = -4;
			// VATParameters::gap_extend = 0;
			// VATParameters::match = 5;
			// VATParameters::penalty = -4;
			po::padding = 8;
			po::hit_cap = 2;
			po::min_identities = 24;
			po::gapped_xdrop = 18;

		}
		else if (po::dna_homology)
		{
			// cout<<"Init DNA homology parameters"<<endl;
			// po::set_option(po::hit_cap, 15u);
			// po::set_option(po::min_identities, 8u);
			// po::set_option(po::gapped_xdrop, 18);
			// po::set_option(po::padding, 8);
			// po::set_option(po::match, 5);
			// po::set_option(po::mismatch, -4);
			// po::set_option(VATParameters::gap_open, 0);
			// po::set_option(VATParameters::gap_extend, 0);
			// po::set_option(po::penalty, -4);
			// po::set_option(po::lowmem, 1u);
			po::gap_open = 0;
			// VATParameters::lowmem = 1;
			po::penalty = -4;
			po::gap_extend = 0;
			po::match = 5;
			po::reward = 5;
			po::penalty = -4;
			po::padding = 8;
			po::hit_cap = 15;
			po::min_identities = 8;
			po::gapped_xdrop = 18;
		}else
		{
			// po::set_option(po::hit_cap, 15u);
			// po::set_option(po::min_identities, 24u);
			// po::set_option(po::gapped_xdrop, 30);
			// po::set_option(po::padding, 8);
			// po::set_option(po::lowmem, 1u);
			po::gap_open = 0;
			// VATParameters::lowmem = 1;
			// VATParameters::penalty = -4;
			po::gap_extend = 0;
			// VATParameters::match = 5;
			// VATParameters::penalty = -4;
			po::padding = 8;
			po::hit_cap = 15u;
			po::min_identities = 24u;
			po::gapped_xdrop = 30;	
		}	
		ScoreMatrix::instance = auto_ptr<ScoreMatrix> (new ScoreMatrix(po::matrix,
				po::gap_open,
				po::gap_extend,
				po::reward,
				po::penalty,
				DNA ()));
		// ScoreMatrix::get().print<DNA>();
	}
	
	verbose_stream << "Gap open penalty = " << po::gap_open << endl;
	verbose_stream << "Gap extension penalty = " << po::gap_extend << endl;

	if(po::seg == "" && po::algn_type == po::blastx)
		po::seg = "yes";
	verbose_stream << "Seg masking = " << (po::seg == "yes") << endl;

	po::have_ssse3 = check_SSSE3();
	if(po::have_ssse3)
		verbose_stream << "SSSE3 enabled." << endl;
	if(po::debug_log) {
		copy_file(log_stream, "/etc/issue");
		copy_file(log_stream, "/proc/cpuinfo");
		copy_file(log_stream, "/proc/meminfo");
	}
	// if (po::chimera)
	// {
	// 	cout<<"Init chimeric alignment parameters"<<endl;
	// 	// po::set_option(po::hit_cap, 15u);
	// 	// po::set_option(po::min_identities, 14u);
	// 	// po::set_option(po::padding, 8);
	// 	// po::set_option(po::match, 1);
	// 	// po::set_option(po::mismatch, -5);
	// 	// po::set_option(po::lowmem, 1u);
	// 	// VATParameters::lowmem = 1u;
	// 	// VATParameters::penalty = -4;
	// 	po::gap_extend = 0;
	// 	po::gap_open = 0;
	// 	po::match = 1;
	// 	po::mismatch = -5;
	// 	po::penalty = -4;
	// 	po::padding = 8;
	// 	po::hit_cap = 15u;
	// 	po::min_identities = 14u;
	// }
	// else if (po::spilce||po::circ)
	// {
	// 	cout<<"Init spliced alignment parameters"<<endl;
	// 	// po::set_option(po::hit_cap, 15u);
	// 	// po::set_option(po::min_identities, 29u);
	// 	// po::set_option(po::padding, 8);
	// 	// po::set_option(po::match, 5);
	// 	// po::set_option(po::mismatch, -2);
	// 	// po::set_option(po::lowmem, 1u);
	// 	// VATParameters::lowmem = 1u;
	// 	// // VATParameters::penalty = -4;
	// 	// // VATParameters::gap_extend = 0;
	// 	po::gap_extend = 1;
	// 	po::gap_open = 1;
	// 	po::match = 5;
	// 	po::mismatch = -2;
	// 	// VATParameters::penalty = -4;
	// 	po::padding = 8;
	// 	po::hit_cap = 15;
	// 	po::min_identities = 29;
	// }
	// else if (po::whole_genome_sequencing)
	// {
	// 	cout<<"Init whole genome sequencing alignment parameters"<<endl;
	// 	// po::set_option(po::hit_cap, 2u);
	// 	// po::set_option(po::min_identities, 28u);
	// 	// po::set_option(po::padding, 8);
	// 	// // po::set_option(po::match, 5);
	// 	// // po::set_option(po::mismatch, -2);
	// 	// po::set_option(po::lowmem, 1u);
	// 	// VATParameters::lowmem = 4;
	// 	// // VATParameters::penalty = -4;
	// 	// // VATParameters::gap_extend = 0;
	// 	po::match = 5;
	// 	po::mismatch = -2;
	// 	// VATParameters::penalty = -4;
	// 	po::padding = 8;
	// 	po::hit_cap = 2;
	// 	po::min_identities = 28;
	// }
	// else if (po::whole_genome)
	// {
	// 	cout<<"Init whole genome alignment parameters"<<endl;
	// 	// po::set_option(po::hit_cap, 2u);
	// 	// po::set_option(po::min_identities, 20u);
	// 	// po::set_option(po::padding, 8);
	// 	// // po::set_option(po::match, 5);
	// 	// // po::set_option(po::mismatch, -2);
	// 	// po::set_option(po::lowmem, 4u);
	// 	// VATParameters::gap_open = 0;
	// 	// VATParameters::lowmem = 4;
	// 	// VATParameters::penalty = -4;
	// 	// VATParameters::gap_extend = 0;
	// 	// VATParameters::match = 5;
	// 	// VATParameters::penalty = -4;
	// 	po::padding = 8;
	// 	po::hit_cap = 2;
	// 	po::min_identities = 24;
	// 	po::gapped_xdrop = 18;

	// }
	// else if (po::dna_homology)
	// {
	// 	cout<<"Init DNA homology parameters"<<endl;
	// 	// po::set_option(po::hit_cap, 15u);
	// 	// po::set_option(po::min_identities, 8u);
	// 	// po::set_option(po::gapped_xdrop, 18);
	// 	// po::set_option(po::padding, 8);
	// 	// po::set_option(po::match, 5);
	// 	// po::set_option(po::mismatch, -4);
	// 	// po::set_option(VATParameters::gap_open, 0);
	// 	// po::set_option(VATParameters::gap_extend, 0);
	// 	// po::set_option(po::penalty, -4);
	// 	// po::set_option(po::lowmem, 1u);
	// 	po::gap_open = 0;
	// 	// VATParameters::lowmem = 1;
	// 	po::penalty = -4;
	// 	po::gap_extend = 0;
	// 	po::match = 5;
	// 	po::penalty = -4;
	// 	po::padding = 8;
	// 	po::hit_cap = 15;
	// 	po::min_identities = 8;
	// 	po::gapped_xdrop = 18;
	// }else
	// {
	// 	// po::set_option(po::hit_cap, 15u);
	// 	// po::set_option(po::min_identities, 24u);
	// 	// po::set_option(po::gapped_xdrop, 30);
	// 	// po::set_option(po::padding, 8);
	// 	// po::set_option(po::lowmem, 1u);
	// 	po::gap_open = 0;
	// 	// VATParameters::lowmem = 1;
	// 	// VATParameters::penalty = -4;
	// 	po::gap_extend = 0;
	// 	// VATParameters::match = 5;
	// 	// VATParameters::penalty = -4;
	// 	po::padding = 8;
	// 	po::hit_cap = 15u;
	// 	po::min_identities = 24u;
	// 	po::gapped_xdrop = 30;	
	// }	
}
void kmer_set()
{
	namespace po = VATParameters;
	if (po::seed_len == 8)
	{
		uint32_t idx = 8 - 8 + 1;
		po::set_option(po::index_mode, idx);
		// cout<<"default = "<<idx<<endl;

	}
	else if (po::seed_len == 9)
	{
		uint32_t idx = 9 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 10)
	{
		uint32_t idx = 10 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 11)
	{
		uint32_t idx = 11 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 12)
	{
		uint32_t idx = 12 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 13)
	{
		uint32_t idx = 13 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 14)
	{
		uint32_t idx = 14 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 16)
	{
		uint32_t idx = 16 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 17)
	{
		uint32_t idx = 17 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 18)
	{
		uint32_t idx = 18 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 19)
	{
		uint32_t idx = 19 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 20)
	{
		uint32_t idx = 20 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 21)
	{
		uint32_t idx = 21 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 22)
	{
		uint32_t idx = 22 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 23)
	{
		uint32_t idx = 23 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 24)
	{
		uint32_t idx = 24 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 25)
	{
		uint32_t idx = 25 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 26)
	{
		uint32_t idx = 26 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 27)
	{
		uint32_t idx = 27 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 28)
	{
		uint32_t idx = 28 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 29)
	{
		uint32_t idx = 29 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 30)
	{
		uint32_t idx = 30 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else if (po::seed_len == 15)
	{
		uint32_t idx = 15 - 8 + 1;
		po::set_option(po::index_mode, idx);
	}
	else
	{
		uint32_t idx = 15 - 8 + 1;
		po::set_option(po::index_mode, idx);
		// cout<<"default = "<<idx<<endl;
	}
}
template<typename _val>
void setup_search_params(pair<size_t,size_t> query_len_bounds, size_t chunk_db_letters)
{
	namespace po = VATParameters;
	if(po::aligner_mode == po::long_model) {
		po::set_option(po::hit_cap, 256u);
	} else if (po::aligner_mode == po::short_model) {
		po::set_option(po::hit_cap, 32u);
	}

	const double b = po::min_bit_score == 0 ? ScoreMatrix::get().bitscore(po::max_evalue, ref_header.letters, query_len_bounds.first) : po::min_bit_score;

	po::set_option(po::min_identities, 18u);
	po::min_ungapped_raw_score = 0;

	po::set_option(po::window, 40u);
	po::set_option(po::hit_band, 5);
	po::min_hit_score = 0;

	// log_stream << "Query len bounds " << query_len_bounds.first << ' ' << query_len_bounds.second << endl;
	// log_stream << "Minimum bit score = " << b << endl;
	// log_stream << "Search parameters " << po::min_ungapped_raw_score << ' ' << po::min_hit_score << ' ' << po::hit_cap << endl;
}

template<>
void setup_search_params<DNA>(pair<size_t,size_t> query_len_bounds, size_t chunk_db_letters)
{

	namespace po = VATParameters;
	if(po::aligner_mode == po::long_model || po::aligner_mode == po::accuracy_model) {
		po::set_option(po::hit_cap, std::max(256u, (unsigned)(chunk_db_letters/8735437)));
	} else if (po::aligner_mode == po::short_model) {
		po::set_option(po::hit_cap, std::max(128u, (unsigned)(chunk_db_letters/17470874)));
	}

	const double b = po::min_bit_score == 0 ? ScoreMatrix::get().bitscore(po::max_evalue, ref_header.letters, query_len_bounds.first) : po::min_bit_score;

	if(query_len_bounds.second <= 40) {
		po::set_option(po::min_identities, 10u);
		po::set_option(po::min_ungapped_raw_score, ScoreMatrix::get().rawscore(std::min(27.0, b)));
	} else {
		po::set_option(po::min_identities, 9u);
		po::set_option(po::min_ungapped_raw_score, ScoreMatrix::get().rawscore(std::min(23.0, b)));
	}

	if(query_len_bounds.second <= 80) {
		const int band = po::read_padding<DNA>(query_len_bounds.second);
		po::set_option(po::window, (unsigned)(query_len_bounds.second + band));
		po::set_option(po::hit_band, band);
		// cout<<"rawscore = "<<score_matrix::get().rawscore(b)<<endl;
		// po::set_option(po::min_hit_score, 11);
		po::set_option(po::min_hit_score, ScoreMatrix::get().rawscore(b));
	} else {
		po::set_option(po::window, 30u);
		po::set_option(po::hit_band, 6);
		// cout<<"rawscore = "<<score_matrix::get().rawscore(std::min(29.0, b))<<endl;
		// po::set_option(po::min_hit_score, 11);
		po::set_option(po::min_hit_score, ScoreMatrix::get().rawscore(std::min(29.0, b)));
	}
	// cout<<VATParameters::hit_cap<<"\t"<<VATParameters::lowmem<<"\t"<<VATParameters::min_identities<<"\t"<<VATParameters::match<<"\t"<<VATParameters::mismatch<<endl;

	// namespace po = VATParameters;
	/*
	if(po::aligner_mode == po::long_model || po::aligner_mode == po::accuracy_model) {
		po::set_option(po::hit_cap, 24u);
	} else if (po::aligner_mode == po::short_model) {
		po::set_option(po::hit_cap, 24u);
	}
*/
	// const double b = 0;
/**
 *         	("band", po::value<int>(&VATParameters::padding)->default_value(0), "band for dynamic programming computation")
       		("max-hits,C", po::value<unsigned>(&VATParameters::hit_cap)->default_value(0), "maximum number of hits to consider for one seed")
       		("id2", po::value<unsigned>(&VATParameters::min_identities)->default_value(25), "minimum number of identities for stage 1 hit")
        	("window,w", po::value<unsigned>(&VATParameters::window)->default_value(0), "window size for local hit search")
        	("xdrop", po::value<int>(&VATParameters::xdrop)->default_value(25), "xdrop for ungapped alignment")
        	("gapped-xdrop,X", po::value<int>(&VATParameters::gapped_xdrop)->default_value(25), "xdrop for gapped alignment in bits")
        	("ungapped-score", po::value<int>(&VATParameters::min_ungapped_raw_score)->default_value(0), "minimum raw alignment score to continue local extension")

	if(query_len_bounds.second <= 40) {
		po::set_option(po::min_identities, 0u);
		po::set_option(po::min_ungapped_raw_score, 0);
	} else {
		po::set_option(po::min_identities, 0u);
		po::set_option(po::min_ungapped_raw_score, 0);
	}

	if(query_len_bounds.second <= 80) {
		const int band = po::read_padding<DNA>(query_len_bounds.second);
		po::set_option(po::window, (unsigned)(query_len_bounds.second + band));
		po::set_option(po::hit_band, band);
		// cout<<"rawscore = "<<score_matrix::get().rawscore(b)<<endl;
		// po::set_option(po::min_hit_score, 11);
		po::set_option(po::min_hit_score, 0);
	} else {
		po::set_option(po::window, 30u);
		po::set_option(po::hit_band, 6);
		// cout<<"rawscore = "<<score_matrix::get().rawscore(std::min(29.0, b))<<endl;
		// po::set_option(po::min_hit_score, 11);
		po::set_option(po::min_hit_score, 0);
	}
	*/
	/*
	if (po::chimera)
	{
		cout<<"Init chimeric alignment parameters"<<endl;
		// po::set_option(po::hit_cap, 15u);
		// po::set_option(po::min_identities, 14u);
		// po::set_option(po::padding, 8);
		// po::set_option(po::match, 1);
		// po::set_option(po::mismatch, -5);
		// po::set_option(po::lowmem, 1u);
		// cout<<VATParameters::hit_cap<<"\t"<<VATParameters::lowmem<<"\t"<<VATParameters::min_identities<<"\t"<<VATParameters::match<<"\t"<<VATParameters::mismatch<<endl;
		// VATParameters::lowmem = 1u;
		// VATParameters::penalty = -4;
		// VATParameters::gap_extend = 0;
		// VATParameters::gap_open = 0;
		// VATParameters::match = 1;
		// VATParameters::mismatch = -5;
		// VATParameters::penalty = -4;
		// VATParameters::padding = 8;
		// VATParameters::hit_cap = 15u;
		// VATParameters::min_identities = 14u;
	}

	if (po::spilce||po::circ)
	{
		cout<<"Init spliced alignment parameters"<<endl;
		// po::set_option(po::hit_cap, 15u);
		// po::set_option(po::min_identities, 29u);
		// po::set_option(po::padding, 8);
		// po::set_option(po::match, 5);
		// po::set_option(po::mismatch, -2);
		// po::set_option(po::lowmem, 1u);
		// VATParameters::lowmem = 1u;
		// // VATParameters::penalty = -4;
		// // VATParameters::gap_extend = 0;
		VATParameters::gap_extend = 1;
		VATParameters::gap_open = 1;
		VATParameters::match = 5;
		VATParameters::mismatch = -2;
		// VATParameters::penalty = -4;
		VATParameters::padding = 8;
		VATParameters::hit_cap = 15;
		VATParameters::min_identities = 29;
	}
	if (po::whole_genome_sequencing)
	{
		cout<<"Init whole genome sequencing alignment parameters"<<endl;
		// po::set_option(po::hit_cap, 2u);
		// po::set_option(po::min_identities, 28u);
		// po::set_option(po::padding, 8);
		// // po::set_option(po::match, 5);
		// // po::set_option(po::mismatch, -2);
		// po::set_option(po::lowmem, 1u);
		// VATParameters::lowmem = 4;
		// // VATParameters::penalty = -4;
		// // VATParameters::gap_extend = 0;
		VATParameters::match = 5;
		VATParameters::mismatch = -2;
		// VATParameters::penalty = -4;
		VATParameters::padding = 8;
		VATParameters::hit_cap = 2;
		VATParameters::min_identities = 28;
	}
	if (po::whole_genome)
	{
		cout<<"Init whole genome alignment parameters"<<endl;
		// po::set_option(po::hit_cap, 2u);
		// po::set_option(po::min_identities, 20u);
		// po::set_option(po::padding, 8);
		// // po::set_option(po::match, 5);
		// // po::set_option(po::mismatch, -2);
		// po::set_option(po::lowmem, 4u);
		// VATParameters::gap_open = 0;
		// VATParameters::lowmem = 4;
		// VATParameters::penalty = -4;
		// VATParameters::gap_extend = 0;
		// VATParameters::match = 5;
		// VATParameters::penalty = -4;
		VATParameters::padding = 8;
		VATParameters::hit_cap = 2;
		VATParameters::min_identities = 24;
		VATParameters::gapped_xdrop = 18;

	}
	if (po::dna_homology)
	{
		cout<<"Init DNA homology parameters"<<endl;
		// po::set_option(po::hit_cap, 15u);
		// po::set_option(po::min_identities, 8u);
		// po::set_option(po::gapped_xdrop, 18);
		// po::set_option(po::padding, 8);
		// po::set_option(po::match, 5);
		// po::set_option(po::mismatch, -4);
		// po::set_option(VATParameters::gap_open, 0);
		// po::set_option(VATParameters::gap_extend, 0);
		// po::set_option(po::penalty, -4);
		// po::set_option(po::lowmem, 1u);
		VATParameters::gap_open = 0;
		// VATParameters::lowmem = 1;
		VATParameters::penalty = -4;
		VATParameters::gap_extend = 0;
		VATParameters::match = 5;
		VATParameters::penalty = -4;
		VATParameters::padding = 8;
		VATParameters::hit_cap = 15;
		VATParameters::min_identities = 8;
		VATParameters::gapped_xdrop = 18;
	}else
	{
		// po::set_option(po::hit_cap, 15u);
		// po::set_option(po::min_identities, 24u);
		// po::set_option(po::gapped_xdrop, 30);
		// po::set_option(po::padding, 8);
		// po::set_option(po::lowmem, 1u);
		VATParameters::gap_open = 0;
		// VATParameters::lowmem = 1;
		// VATParameters::penalty = -4;
		VATParameters::gap_extend = 0;
		// VATParameters::match = 5;
		// VATParameters::penalty = -4;
		VATParameters::padding = 8;
		VATParameters::hit_cap = 15u;
		VATParameters::min_identities = 24u;
		VATParameters::gapped_xdrop = 30;	
	}	
	*/			
}


template<>
void setup_search_params<Protein>(pair<size_t,size_t> query_len_bounds, size_t chunk_db_letters)
{
	namespace po = VATParameters;
	if(po::aligner_mode == po::long_model || po::aligner_mode == po::accuracy_model) {
		po::set_option(po::hit_cap, std::max(256u, (unsigned)(chunk_db_letters/8735437)));
	} else if (po::aligner_mode == po::short_model) {
		po::set_option(po::hit_cap, std::max(128u, (unsigned)(chunk_db_letters/17470874)));
	}

	const double b = po::min_bit_score == 0 ? ScoreMatrix::get().bitscore(po::max_evalue, ref_header.letters, query_len_bounds.first) : po::min_bit_score;

	if(query_len_bounds.second <= 40) {
		po::set_option(po::min_identities, 10u);
		po::set_option(po::min_ungapped_raw_score, ScoreMatrix::get().rawscore(std::min(27.0, b)));
	} else {
		po::set_option(po::min_identities, 9u);
		po::set_option(po::min_ungapped_raw_score, ScoreMatrix::get().rawscore(std::min(23.0, b)));
	}

	if(query_len_bounds.second <= 80) {
		const int band = po::read_padding<Protein>(query_len_bounds.second);
		po::set_option(po::window, (unsigned)(query_len_bounds.second + band));
		po::set_option(po::hit_band, band);
		po::set_option(po::min_hit_score, ScoreMatrix::get().rawscore(b));
	} else {
		po::set_option(po::window, 40u);
		po::set_option(po::hit_band, 5);
		po::set_option(po::min_hit_score, ScoreMatrix::get().rawscore(std::min(29.0, b)));
	}
	// log_stream << "Query len bounds " << query_len_bounds.first << ' ' << query_len_bounds.second << endl;
	// log_stream << "Search parameters " << po::min_ungapped_raw_score << ' ' << po::min_hit_score << ' ' << po::hit_cap << endl;
}

#endif /* SETUP_H_ */
