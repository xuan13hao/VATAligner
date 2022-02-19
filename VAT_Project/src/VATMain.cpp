

#include <iostream>
#include <iterator>
#include <boost/program_options.hpp>
#include "basic/VATParameters.h"
#include "util/log_stream.h"
#include "data/reference.h"
#include "run/CreateDB.h"
#include "run/AlignmentFlow.h"
#include "util/complexity_filter.h"
#include "basic/ContextSet.h"
#include "output/view.h"


using std::cout;
using std::cerr;
using std::endl;

int main(int ac, const char* av[])
{

	namespace po = boost::program_options;

	try {

		string command;

        po::options_description general("General options");
        general.add_options()
            ("help,h", "produce help message")
            ("threads,p", po::value<uint32_t>(&VATParameters::threads_)->default_value(0), "number of cpu threads")
            ("db,d", po::value<string>(&VATParameters::database), "database file")
            ("daa,a", po::value<string>(&VATParameters::daa_file), "VAT alignment archive (DAA) file")
            ("verbose,v", "enable verbose out")
            ("log", "enable debug log");

        po::options_description makedb("Makedb options");
        makedb.add_options()
        	("in", po::value<string>(&VATParameters::input_ref_file), "input reference file in FASTA format")
        	("block-size,b", po::value<double>(&VATParameters::chunk_size), "sequence block size in billions of letters (default=2)")
#ifdef EXTRA
        	("dbtype", po::value<string>(&program_options::db_type), "database type (nucl/prot)")
#endif
        	;
        	//("kegg-map", po::value<string>(&program_options::kegg_file), "KEGG mapping file");

        po::options_description aligner("Aligner options");
        aligner.add_options()
			("query,q", po::value<string>(&VATParameters::query_file), "input query file")
			("max-target-seqs,k", po::value<uint64_t>(&VATParameters::max_alignments)->default_value(25), "maximum number of target sequences to report alignments for")
			("top", po::value<double>(&VATParameters::toppercent)->default_value(100), "report alignments within this percentage range of top alignment score (overrides --max-target-seqs)")
        	("compress", po::value<unsigned>(&VATParameters::compression)->default_value(0), "compression for output files (0=none, 1=gzip)")
			("evalue,e", po::value<double>(&VATParameters::max_evalue)->default_value(0.001), "maximum e-value to report alignments")
        	("min-score", po::value<double>(&VATParameters::min_bit_score)->default_value(0), "minimum bit score to report alignments (overrides e-value setting)")
        	("id", po::value<double>(&VATParameters::min_id)->default_value(0), "minimum identity% to report an alignment")
        	("sensitive", "enable sensitive mode (default: fast)")
        	("index-chunks,c", po::value<unsigned>(&VATParameters::lowmem)->default_value(4), "number of chunks for index processing")
        	("tmpdir,t", po::value<string>(&VATParameters::tmpdir)->default_value("/dev/shm"), "directory for temporary files")
        	("gapopen", po::value<int>(&VATParameters::gap_open)->default_value(-1), "gap open penalty, -1=default (11 for protein)")
        	("gapextend", po::value<int>(&VATParameters::gap_extend)->default_value(-1), "gap extension penalty, -1=default (1 for protein)")
#ifdef EXTRA
        	("reward", po::value<int>(&program_options::reward)->default_value(2), "match reward score (blastn only)")
        	("penalty", po::value<int>(&program_options::penalty)->default_value(-3), "mismatch penalty score (blastn only)")
#endif
        	("matrix", po::value<string>(&VATParameters::matrix)->default_value("blosum62"), "score matrix for protein alignment")
        	("seg", po::value<string>(&VATParameters::seg), "enable SEG masking of queries (yes/no)");
			//("salltitles", "print all subject titles into the blast tabular format");
        	//("very-sensitive", "enable very sensitive mode (default: fast)");

        po::options_description advanced("Advanced options (0=auto)");
        advanced.add_options()
			("seed-freq", po::value<double>(&VATParameters::max_seed_freq)->default_value(-15), "maximum seed frequency")
			("run-len,l", po::value<unsigned>(&VATParameters::run_len)->default_value(0), "mask runs between stop codons shorter than this length")
       		("max-hits,C", po::value<unsigned>(&VATParameters::hit_cap)->default_value(0), "maximum number of hits to consider for one seed")
       		("id2", po::value<unsigned>(&VATParameters::min_identities)->default_value(0), "minimum number of identities for stage 1 hit")
        	("window,w", po::value<unsigned>(&VATParameters::window)->default_value(0), "window size for local hit search")
        	("xdrop", po::value<int>(&VATParameters::xdrop)->default_value(20), "xdrop for ungapped alignment")
        	("gapped-xdrop,X", po::value<int>(&VATParameters::gapped_xdrop)->default_value(20), "xdrop for gapped alignment in bits")
        	("ungapped-score", po::value<int>(&VATParameters::min_ungapped_raw_score)->default_value(0), "minimum raw alignment score to continue local extension")
        	("hit-band", po::value<int>(&VATParameters::hit_band)->default_value(0), "band for hit verification")
        	("hit-score", po::value<int>(&VATParameters::min_hit_score)->default_value(0), "minimum score to keep a tentative alignment")
        	("band", po::value<int>(&VATParameters::padding)->default_value(0), "band for dynamic programming computation")
        	("shapes,s", po::value<unsigned>(&VATParameters::shapes)->default_value(0), "number of seed shapes (0 = all available)")
        	("index-mode", po::value<unsigned>(&VATParameters::index_mode)->default_value(0), "index mode (1=4x12, 2=16x9)")
        	("fetch-size", po::value<unsigned>(&VATParameters::fetch_size)->default_value(4096), "trace point fetch size")
        	("single-domain", "Discard secondary domains within one target sequence")
        	("no-traceback,r", "disable alignment traceback")
        	("dbsize", po::value<size_t>(&VATParameters::db_size)->default_value(0), "effective database size (in letters)");
        	//("compress-temp", po::value<unsigned>(&program_options::compress_temp)->default_value(0), "compression for temporary output files (0=none, 1=gzip)");

        po::options_description view_options("View options");
        view_options.add_options()
			("out,o", po::value<string>(&VATParameters::output_file), "output file")
			("outfmt,f", po::value<string>(&VATParameters::output_format)->default_value("tab"), "output format (tab/sam)")
			("forwardonly", "only show alignments of forward strand");

        po::options_description hidden("Hidden options");
        hidden.add_options()
        	("command", po::value<string>(&command))
#ifdef EXTRA
        	("match1", po::value<string>(&program_options::match_file1))
        	("match2", po::value<string>(&program_options::match_file2))
        	("tab", "tabular format")
#endif
        	;

        po::options_description cmd_line_options("Command line options");
        cmd_line_options.add(general).add(hidden).add(makedb).add(aligner).add(advanced).add(view_options);

        po::positional_options_description positional;
        positional.add("command", -1);

        po::variables_map vm;
        po::store(po::command_line_parser(ac, av).options(cmd_line_options).positional(positional).run(), vm);
        po::notify(vm);

        if(vm.count("sensitive"))
        	VATParameters::aligner_mode = VATParameters::sensitive;
        //else if(vm.count("very-sensitive"))
        //	program_options::aligner_mode = program_options::very_sensitive;
        else
        	VATParameters::aligner_mode = VATParameters::fast;
        VATParameters::alignment_traceback = (vm.count("no-traceback") == 0);
        VATParameters::long_mode = vm.count("long") > 0;
        VATParameters::verbose = vm.count("verbose") > 0;
        VATParameters::debug_log = vm.count("log") > 0;
        VATParameters::salltitles = vm.count("salltitles") > 0;
        VATParameters::forwardonly = vm.count("forwardonly") > 0;
        VATParameters::single_domain = vm.count("single-domain") > 0;

        setup(command, ac, av);

        if (vm.count("help")) {
        	cout << endl << "Syntax:" << endl;
        	cout << "  diamond COMMAND [OPTIONS]" << endl << endl;
        	cout << "Commands:" << endl;
        	cout << "  makedb\tBuild diamond database from a FASTA file" << endl;
        	cout << "  blastp\tAlign amino acid query sequences against a protein reference database" << endl;
        	cout << "  blastx\tAlign DNA query sequences against a protein reference database" << endl;
        	cout << "  view\tView DIAMOND alignment archive (DAA) formatted file" << endl;
        	cout << endl;
        	cout << general << endl << makedb << endl << aligner << endl << advanced << endl << view_options << endl;
        } else if (VATParameters::command == VATParameters::makedb && vm.count("in") && vm.count("db")) 
		{
        	if(vm.count("block-size") == 0)
        		VATParameters::chunk_size = 2;
        		// make_db(Amino_acid());
				// make_db(DNA());
				// make_db(Protein());


        } else if ((VATParameters::command == VATParameters::blastp
        		|| VATParameters::command == VATParameters::blastx
				|| VATParameters::command == VATParameters::blastn
#ifdef EXTRA
        		|| program_options::command == program_options::blastn
#endif
        		)
        		&& vm.count("query") && vm.count("db") && vm.count("daa")) 
		{
        	if(vm.count("block-size") > 0) {
        		cerr << "Warning: --block-size option should be set for the makedb command." << endl;
        	} else
        		VATParameters::chunk_size = 0;
        	if(VATParameters::command == VATParameters::blastp)
				master_thread<Protein>();
        		// master_thread<DNA>();
        	// else if(program_options::command == program_options::blastp)
        	// 	master_thread<Protein>();
        } else if(VATParameters::command == VATParameters::view && vm.count("daa") > 0)
        	view();
		#ifdef EXTRA
        else if (command == "stat" && vm.count("match1"))
        	if(program_options::db_type == "nucl")
        		blast_stat<Nucleotide>(vm.count("tab") > 0);
        	else
        		blast_stat<Amino_acid>(vm.count("tab") > 0);
        //else if (command == "comp" && vm.count("query") && vm.count("match1") && vm.count("match2"))
        	//compare();
        //else if (command == "random" && vm.count("query") && vm.count("out"))
        	//random_seq();
		#endif
		#ifdef EXTRA
        else if (command == "test")
        	test_io();
		#endif
        else
        	cout << "Insufficient arguments. Use diamond -h for help.\n";
	}
	catch(std::bad_alloc &e) {
		cerr << "Failed to allocate sufficient memory. Please refer to the readme for instructions on memory usage." << endl;
		log_stream << "Error: " << e.what() << endl;
	} catch(exception& e) {
        cerr << "Error: " << e.what() << endl;
        log_stream << "Error: " << e.what() << endl;
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        return 1;
    }

    return 0;
}
