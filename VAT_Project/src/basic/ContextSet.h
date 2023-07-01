
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
			po::gap_open = 3;
		if(po::gap_extend == -1)
			po::gap_extend = 2;
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
	po::min_ungapped_raw_score = ScoreMatrix::get().rawscore(std::min(po::min_ungapped_raw_score == 0 ? 19.0 : po::min_ungapped_raw_score, b));

	po::set_option(po::window, 40u);
	po::set_option(po::hit_band, 5);
	po::min_hit_score = ScoreMatrix::get().rawscore(std::min(po::min_hit_score == 0 ? 19.0 : po::min_hit_score, b));

	log_stream << "Query len bounds " << query_len_bounds.first << ' ' << query_len_bounds.second << endl;
	log_stream << "Minimum bit score = " << b << endl;
	log_stream << "Search parameters " << po::min_ungapped_raw_score << ' ' << po::min_hit_score << ' ' << po::hit_cap << endl;
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
	log_stream << "Query len bounds " << query_len_bounds.first << ' ' << query_len_bounds.second << endl;
	log_stream << "Search parameters " << po::min_ungapped_raw_score << ' ' << po::min_hit_score << ' ' << po::hit_cap << endl;
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
	log_stream << "Query len bounds " << query_len_bounds.first << ' ' << query_len_bounds.second << endl;
	log_stream << "Search parameters " << po::min_ungapped_raw_score << ' ' << po::min_hit_score << ' ' << po::hit_cap << endl;
}

#endif /* SETUP_H_ */
