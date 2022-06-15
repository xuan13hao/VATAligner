

#include <iostream>
#include <iterator>
#include <boost/program_options.hpp>
#include "basic/VATParameters.h"
#include "tools/TimerTools.h"
#include "data/Reference.h"
#include "model/RunModel.h"
#include "tools/complexity_filter.h"
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
            ("vaa,a", po::value<string>(&VATParameters::daa_file), "VAT alignment archive (vaa) file")
            ("verbose,v", "enable verbose out")
            ("log", "enable debug log")
        	("dbtype", po::value<string>(&VATParameters::db_type), "database type (nucl/prot)");

        po::options_description makedb("Makedb options");
        makedb.add_options()
        	("in", po::value<string>(&VATParameters::input_ref_file), "input reference file in FASTA format")
        	("block-size,b", po::value<double>(&VATParameters::chunk_size), "sequence block size in billions of letters (default=2)")
        	;

        po::options_description aligner("Aligner options");
        aligner.add_options()
			("query,q", po::value<string>(&VATParameters::query_file), "input query file")
			("max-target-seqs,k", po::value<uint64_t>(&VATParameters::max_alignments)->default_value(25), "maximum number of target sequences to report alignments for")
			("top", po::value<double>(&VATParameters::toppercent)->default_value(100), "report alignments within this percentage range of top alignment score (overrides --max-target-seqs)")
        	("compress", po::value<unsigned>(&VATParameters::compression)->default_value(0), "compression for output files (0=none, 1=gzip)")
			("evalue,e", po::value<double>(&VATParameters::max_evalue)->default_value(0.001), "maximum e-value to report alignments")
        	("min-score", po::value<double>(&VATParameters::min_bit_score)->default_value(0), "minimum bit score to report alignments (overrides e-value setting)")
        	("id", po::value<double>(&VATParameters::min_id)->default_value(30), "minimum identity% to report an alignment")
        	("sensitive", "enable sensitive mode (default: fast)")
        	("index-chunks,c", po::value<unsigned>(&VATParameters::lowmem)->default_value(4), "number of chunks for index processing")
        	("tmpdir,t", po::value<string>(&VATParameters::tmpdir)->default_value("/dev/shm"), "directory for temporary files")
        	("gapopen", po::value<int>(&VATParameters::gap_open)->default_value(-1), "gap open penalty, -1=default (11 for protein)")
        	("gapextend", po::value<int>(&VATParameters::gap_extend)->default_value(-1), "gap extension penalty, -1=default (1 for protein)")
        	("reward", po::value<int>(&VATParameters::reward)->default_value(2), "match reward score (blastn only)")
        	("penalty", po::value<int>(&VATParameters::penalty)->default_value(-3), "mismatch penalty score (blastn only)")

        	("matrix", po::value<string>(&VATParameters::matrix)->default_value("blosum62"), "score matrix for protein alignment")
        	("seg", po::value<string>(&VATParameters::seg), "enable SEG masking of queries (yes/no)");

        po::options_description advanced("Advanced options (0=auto)");
        advanced.add_options()
			("seed-freq", po::value<double>(&VATParameters::max_seed_freq)->default_value(-15), "maximum seed frequency")
			("run-len,l", po::value<unsigned>(&VATParameters::run_len)->default_value(0), "mask runs between stop codons shorter than this length")
       		("max-hits,C", po::value<unsigned>(&VATParameters::hit_cap)->default_value(0), "maximum number of hits to consider for one seed")
       		("id2", po::value<unsigned>(&VATParameters::min_identities)->default_value(30), "minimum number of identities for stage 1 hit")
        	("window,w", po::value<unsigned>(&VATParameters::window)->default_value(0), "window size for local hit search")
        	("xdrop", po::value<int>(&VATParameters::xdrop)->default_value(20), "xdrop for ungapped alignment")
        	("gapped-xdrop,X", po::value<int>(&VATParameters::gapped_xdrop)->default_value(10), "xdrop for gapped alignment in bits")
        	("ungapped-score", po::value<int>(&VATParameters::min_ungapped_raw_score)->default_value(0), "minimum raw alignment score to continue local extension")
        	("hit-band", po::value<int>(&VATParameters::hit_band)->default_value(0), "band for hit verification")
        	("hit-score", po::value<int>(&VATParameters::min_hit_score)->default_value(0), "minimum score to keep a tentative alignment")
        	("band", po::value<int>(&VATParameters::padding)->default_value(0), "band for dynamic programming computation")
        	("shapes,s", po::value<unsigned>(&VATParameters::shapes)->default_value(0), "number of seed shapes (0 = all available)")
        	("index-mode", po::value<unsigned>(&VATParameters::index_mode)->default_value(0), "index mode")//future interface
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

        if (vm.count("help")) 
		{
        	cout << endl << "Syntax:" << endl;
        	cout << "  VAT COMMAND [OPTIONS]" << endl << endl;
        	cout << "Commands:" << endl;
        	cout << "  makevatdb\tBuild VAT database from a FASTA file" << endl;
        	cout << "  protein\tAlign amino acid query sequences against a protein reference database" << endl;
        	cout << "  dna\tAlign DNA query sequences against a DNA reference database" << endl;
        	cout << "  view\tView VAT alignment archive (vaa) formatted file" << endl;
        	cout << endl;
        	cout << general << endl << makedb << endl << aligner << endl << advanced << endl << view_options << endl;
        } 
		else if (VATParameters::command == VATParameters::makevatdb && vm.count("in") && vm.count("db")) 
		{
        	if(vm.count("block-size") == 0)
			{
				VATParameters::chunk_size = 2;
				// RunModel::CreateDNADB();

				if (vm.count("dbtype")&&VATParameters::db_type == "nucl")
				{
					RunModel::CreateDNADB();
				}
				else if (vm.count("dbtype")&&VATParameters::db_type == "prot")
				{
					RunModel::CreateProteinDB();
				}else
				{
					cerr << "Failed to get databasetype Please refer to the readme for instructions." << endl;
				}
				
			}

        } else if ((VATParameters::command == VATParameters::protein
        		//|| VATParameters::command == VATParameters::blastx
				|| VATParameters::command == VATParameters::dna)
        		&& vm.count("query") && vm.count("db") && vm.count("vaa")) 
		{
        	if(vm.count("block-size") > 0) 
			{
        		cerr << "Warning: --block-size option should be set for the makevatdb command." << endl;
        	} else
			{
				VATParameters::chunk_size = 0;
				if(VATParameters::command == VATParameters::protein)
				{
					cout<<"protein"<<endl;
					RunModel::ProteinAlign();
				}else if (VATParameters::command == VATParameters::dna)
				{
					cout<<"dna"<<endl;
					RunModel::DNAAlign();
				}else
				{
					cerr << "Failed to get alignment type. Please refer to the readme for instructions." << endl;
				}
			}
        } else if(VATParameters::command == VATParameters::view && vm.count("vaa") > 0)
        	view();
        else
        	cout << "Insufficient arguments. Use VAT -h for help.\n";
	}
	catch(std::bad_alloc &e) 
	{
		cerr << "Failed to allocate sufficient memory. Please refer to the readme for instructions on memory usage." << endl;
		log_stream << "Error: " << e.what() << endl;
	} catch(exception& e) 
	{
        cerr << "Error: " << e.what() << endl;
        log_stream << "Error: " << e.what() << endl;
        return 1;
    }
    catch(...) 
	{
        cerr << "Exception of unknown type!\n";
        return 1;
    }

    return 0;
}
