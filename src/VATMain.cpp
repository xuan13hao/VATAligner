

#include <iostream>
#include <iterator>
#include <boost/program_options.hpp>
#include "Commons/VATParameters.h"
#include "tools/TimerTools.h"
#include "Database/Reference.h"
#include "model/RunModel.h"
#include "tools/complexity_filter.h"
#include "Commons/ContextSet.h"
#include "Output/view.h"


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
            ("threads,p", po::value<uint32_t>(&VATParameters::threads_)->default_value(4), "number of cpu threads")
            ("db,d", po::value<string>(&VATParameters::database), "database file")
            ("vaa,a", po::value<string>(&VATParameters::daa_file), "VAT alignment archive (vatr) file")
        	("dbtype", po::value<string>(&VATParameters::db_type), "database type (nucl/prot)");

        po::options_description makedb("Makedb options");
        makedb.add_options()
        	("in,i", po::value<string>(&VATParameters::input_ref_file), "input reference file in FASTA format")
        	// ("block-size,b", po::value<double>(&VATParameters::chunk_size)->default_value(4), "sequence block size in billions of letters (default=4)")
        	;

        po::options_description aligner("Aligner options");
        aligner.add_options()
			("query,q", po::value<string>(&VATParameters::query_file), "input query file")
			("maxtarget_seqs,k", po::value<uint64_t>(&VATParameters::max_alignments)->default_value(25), "maximum number of target sequences to report alignments for")
			("top", po::value<double>(&VATParameters::toppercent)->default_value(98), "report alignments within this percentage range of top alignment score (overrides --maxtarget-seqs)")
			("evalue,e", po::value<double>(&VATParameters::max_evalue)->default_value(0.001), "maximum e-value to report alignments")
        	("min_score", po::value<double>(&VATParameters::min_bit_score)->default_value(0), "minimum bit score to report alignments (overrides e-value setting)")
        	("report_id", po::value<double>(&VATParameters::min_id)->default_value(0), "minimum identity% to report an alignment")
        	("chunks", po::value<unsigned>(&VATParameters::lowmem)->default_value(1), "number of chunks for index processing")
        	("tmpdir,t", po::value<string>(&VATParameters::tmpdir)->default_value("/dev/shm"), "directory for temporary files")
        	("gapopen", po::value<int>(&VATParameters::gap_open)->default_value(-1), "gap open penalty, -1=default (11 for protein)")
        	("gapextend", po::value<int>(&VATParameters::gap_extend)->default_value(-1), "gap extension penalty, -1=default (1 for protein)")
        	("reward", po::value<int>(&VATParameters::reward)->default_value(2), "match reward score (blastn only)")
			("seed_len,S", po::value<int>(&VATParameters::seed_len)->default_value(15), "Seed length default(15) for DNA, Seed length default(8) for Protein")
        	("penalty", po::value<int>(&VATParameters::penalty)->default_value(-3), "mismatch penalty score (blastn only)")
			("match", po::value<int>(&VATParameters::match)->default_value(5), "match score (5)")
			("mismatch", po::value<int>(&VATParameters::mismatch)->default_value(-4), "mismatch score (-4)")
			("simd_sort", "Enable Double Index Sorting on AVX2")
			("chimera", "Enable Chimera alignment")
			("circ", "Enable Circ alignment")
			("wga", "Enable Whole-genome alignment")
			("wgs", "Enable Whole-genome sequencing")
			("splice", "Enable Splice alignments ")
			("dnah", "Enable DNA homology ")
			("avx2", "Enable AVX2 hamming distance ")
			("SEN", "Enable sensitive protein alignment")
        	("matrix", po::value<string>(&VATParameters::matrix)->default_value("blosum62"), "score matrix for protein alignment")
			;
        po::options_description advanced("Advanced options (0=auto)");
        advanced.add_options()
			// ("seed_freq", po::value<double>(&VATParameters::max_seed_freq)->default_value(-20), "maximum seed frequency")
       		("max_seeds,M", po::value<unsigned>(&VATParameters::hit_cap)->default_value(0), "maximum number of hits to consider for one seed")
        	("window,W", po::value<unsigned>(&VATParameters::window)->default_value(0), "window size for local hit search")
			("minimizer,w", po::value<unsigned>(&VATParameters::mini_mizer)->default_value(10), "window size for minimizer")
        	("xdrop", po::value<int>(&VATParameters::xdrop)->default_value(18), "xdrop for ungapped alignment")
        	("gapped_xdrop,X", po::value<int>(&VATParameters::gapped_xdrop)->default_value(18), "xdrop for gapped alignment in bits")
        	("ungapped_score", po::value<int>(&VATParameters::min_ungapped_raw_score)->default_value(0), "minimum raw alignment score to continue local extension")
        	// ("hit-band", po::value<int>(&VATParameters::hit_band)->default_value(0), "band for hit verification")
        	("pre_score", po::value<int>(&VATParameters::min_hit_score)->default_value(0), "minimum score to keep a pre-alignment")
        	("band", po::value<int>(&VATParameters::padding)->default_value(8), "band for dynamic programming computation")
        	("num_shapes,N", po::value<unsigned>(&VATParameters::shapes)->default_value(0), "number of seed shapes (0 = all available)")
			("spaced", po::value<string>(&VATParameters::spaced_seed)->default_value("null"), "Spaced seed")
			("ra", po::value<string>(&VATParameters::reduced_alphabet)->default_value("null"), "Reduced alphabet patterns (dssp.5, murphy.5, dssp.10, murphy.10, MMSEQS12,and td.10)")
        	// ("mode_shape", po::value<unsigned>(&VATParameters::index_mode)->default_value(0), "index mode")//future interface
			("for_only", "only forward strand alignment")
			;
	


        po::options_description view_options("View options");
        view_options.add_options()
			("out,o", po::value<string>(&VATParameters::output_file), "output file")
			("outfmt,f", po::value<string>(&VATParameters::output_format)->default_value("tab"), "output format (tab/sam/paf)")
			// ("forwardonly", "only show alignments of forward strand")
			;

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
        if(vm.count("SEN"))
        	VATParameters::aligner_mode = VATParameters::accuracy_model;
		VATParameters::lowmem = 1;
		VATParameters::chunk_size = 4;
        VATParameters::alignment_traceback = (vm.count("no-traceback") == 0);
        // VATParameters::long_mode = vm.count("long") > 0;
        // VATParameters::verbose = vm.count("verbose") > 0;
        VATParameters::debug_log = vm.count("log") > 0;
        VATParameters::salltitles = vm.count("salltitles") > 0;
        VATParameters::forwardonly = vm.count("forwardonly") > 0;
		VATParameters::forward_only = vm.count("for_only") > 0;
		VATParameters::chimera = vm.count("chimera") > 0;
		VATParameters::whole_genome_sequencing = vm.count("wgs") > 0;
		VATParameters::dna_homology = vm.count("dnah") > 0;
		VATParameters::whole_genome = vm.count("wga") > 0;
		VATParameters::circ = vm.count("circ") > 0;
		VATParameters::spilce = vm.count("splice") > 0;
        // VATParameters::single_domain = vm.count("single-domain") > 0;
		VATParameters::enable_avx2 = vm.count("avx2") > 0;

        setup(command, ac, av);

        if (vm.count("help")) 
		{
        	cout << endl << "Syntax:" << endl;
        	cout << "  VAT COMMAND [OPTIONS]" << endl << endl;
        	cout << "Commands:" << endl;
        	cout << "makevatdb\tBuild VAT database from a FASTA file" << endl;
        	cout << "protein\tAlign protein query sequences against a protein reference database" << endl;
        	cout << "dna\tAlign DNA query sequences against a DNA reference database" << endl;
			// cout <<"blastx\tAlign DNA query sequences against a protein reference database" << endl;
        	cout << "  view\tView VAT alignment archive (vaa) formatted file" << endl;
        	cout << endl;
        	cout << general << endl << makedb << endl << aligner << endl << advanced << endl << view_options << endl;
        } 
		else if (VATParameters::algn_type == VATParameters::makevatdb && vm.count("in") && vm.count("db")) 
		{
        	// if(vm.count("block-size") == 0)
			// {
				if (vm.count("dbtype")&&VATParameters::db_type == "nucl")
				{
					VATParameters::chunk_size = 8;
					RunModel::CreateDNADB();
				}
				else if (vm.count("dbtype")&&VATParameters::db_type == "prot")
				{
					VATParameters::chunk_size = 8;
					RunModel::CreateProteinDB();
				}else
				{
					cerr << "Failed to get databasetype. Please refer to the readme for instructions." << endl;
				}
				
			// }

        } else if ((VATParameters::algn_type == VATParameters::protein
        		|| VATParameters::algn_type == VATParameters::blastx
				|| VATParameters::algn_type == VATParameters::dna)
        		&& vm.count("query") && vm.count("db") && vm.count("vaa")) 
		{
        	// if(vm.count("block-size") > 0) 
			// {
        	// 	cerr << "Warning: --block-size option should be set for the makevatdb command." << endl;
        	// } else
			// {
				if(VATParameters::algn_type == VATParameters::protein)
				{
					// VATParameters::chunk_size = 8;
					RunModel::ProteinAlign();
				}else if (VATParameters::algn_type == VATParameters::dna)
				{
					// VATParameters::chunk_size = 8;
					RunModel::DNAAlign();
				}else if (VATParameters::algn_type == VATParameters::blastx)
				{
					RunModel::BlastxAlign();
				}else
				{
					cerr << "Failed to get alignment type. Please refer to the readme for instructions." << endl;
				}
			// }
        } else if(VATParameters::algn_type == VATParameters::view && vm.count("vaa") > 0)
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
