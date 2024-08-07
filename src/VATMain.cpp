

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
            ("threads,p", po::value<uint32_t>(&VATParameters::threads_)->default_value(1), "number of cpu threads")
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
			("top", po::value<double>(&VATParameters::toppercent)->default_value(98), "report alignments within this percentage range of top alignment score (overrides --max-target-seqs)")
        	// ("compress", po::value<unsigned>(&VATParameters::compression)->default_value(0), "compression for output files (0=none, 1=gzip)")
			("evalue,e", po::value<double>(&VATParameters::max_evalue)->default_value(0.001), "maximum e-value to report alignments")
        	("min_score", po::value<double>(&VATParameters::min_bit_score)->default_value(0), "minimum bit score to report alignments (overrides e-value setting)")
        	("report_id", po::value<double>(&VATParameters::min_id)->default_value(0), "minimum identity% to report an alignment")
        	// ("long-read", "enable long-read mode (default: short)")
			// ("accuracy", "enable accuracy mode")
        	// ("index-chunks,c", po::value<unsigned>(&VATParameters::lowmem)->default_value(1), "number of chunks for index processing")
        	("tmpdir,t", po::value<string>(&VATParameters::tmpdir)->default_value("/dev/shm"), "directory for temporary files")
        	("gapopen", po::value<int>(&VATParameters::gap_open)->default_value(-1), "gap open penalty, -1=default (11 for protein)")
        	("gapextend", po::value<int>(&VATParameters::gap_extend)->default_value(-1), "gap extension penalty, -1=default (1 for protein)")
        	("reward", po::value<int>(&VATParameters::reward)->default_value(2), "match reward score (blastn only)")
			("seed_len,S", po::value<int>(&VATParameters::seed_len)->default_value(15), "Seed length default(15) for DNA, Seed length default(8) for Protein")
        	("penalty", po::value<int>(&VATParameters::penalty)->default_value(-3), "mismatch penalty score (blastn only)")
			("match", po::value<int>(&VATParameters::match)->default_value(5), "match score (5)")
			("mismatch", po::value<int>(&VATParameters::mismatch)->default_value(-4), "mismatch score (-4)")
			// ("whole-genome", po::value<bool>(&VATParameters::whole_genome)->default_value(0), "whole genome alignment (0)")
			("simd_sort", "Double-index based on SIMD")
			("chimera", "Chimera alignment")
			("circ", "Circ alignment")
			("wga", "Whole-genome alignment")
			("wgs", "Whole-genome sequencing")
			("splice", "Splice alignments ")
			("dnah", "DNA homology ")
			("spaced", po::value<string>(&VATParameters::spaced_seed)->default_value("null"), "Spaced seed")
        	("matrix", po::value<string>(&VATParameters::matrix)->default_value("blosum62"), "score matrix for protein alignment")
        	// ("seg", po::value<string>(&VATParameters::seg), "enable SEG masking of queries (yes/no)")
			;
//1111111111111111
        po::options_description advanced("Advanced options (0=auto)");
        advanced.add_options()
			// ("seed_freq", po::value<double>(&VATParameters::max_seed_freq)->default_value(-20), "maximum seed frequency")
			// ("run-len,l", po::value<unsigned>(&VATParameters::run_len)->default_value(0), "mask runs between stop codons shorter than this length")
       		("max_hits,C", po::value<unsigned>(&VATParameters::hit_cap)->default_value(0), "maximum number of hits to consider for one seed")
       		("pre_filter", po::value<unsigned>(&VATParameters::min_identities)->default_value(0), "minimum number of identities for pre-filter hit")
        	// ("window,w", po::value<unsigned>(&VATParameters::window)->default_value(0), "window size for local hit search")
        	("xdrop", po::value<int>(&VATParameters::xdrop)->default_value(18), "xdrop for ungapped alignment")
        	("gapped_xdrop,X", po::value<int>(&VATParameters::gapped_xdrop)->default_value(18), "xdrop for gapped alignment in bits")
        	("ungapped_score", po::value<int>(&VATParameters::min_ungapped_raw_score)->default_value(0), "minimum raw alignment score to continue local extension")
        	// ("hit-band", po::value<int>(&VATParameters::hit_band)->default_value(0), "band for hit verification")
        	("hit_score", po::value<int>(&VATParameters::min_hit_score)->default_value(0), "minimum score to keep a tentative alignment")
        	("band", po::value<int>(&VATParameters::padding)->default_value(8), "band for dynamic programming computation")
        	// ("shapes,s", po::value<unsigned>(&VATParameters::shapes)->default_value(0), "number of seed shapes (0 = all available)")
        	// ("index-mode", po::value<unsigned>(&VATParameters::index_mode)->default_value(0), "index mode")//future interface
        	// ("fetch-size", po::value<unsigned>(&VATParameters::fetch_size)->default_value(4096), "trace point fetch size")
        	// ("single-domain", "Discard secondary domains within one target sequence")
        	// ("no-traceback,r", "disable alignment traceback")
			("for_only", "only forward strand alignment")
        	// ("dbsize", po::value<size_t>(&VATParameters::db_size)->default_value(0), "effective database size (in letters)")
			;
	
        	//("compress-temp", po::value<unsigned>(&program_options::compress_temp)->default_value(0), "compression for temporary output files (0=none, 1=gzip)");

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

        // if(vm.count("long-read"))
		// {
		// 	VATParameters::aligner_mode = VATParameters::long_model;
		// }
		// else if (vm.count("accuracy"))
		// {
		// 	VATParameters::aligner_mode = VATParameters::accuracy_model;
		// }
        // else
        // 	VATParameters::aligner_mode = VATParameters::short_model;
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
        VATParameters::single_domain = vm.count("single-domain") > 0;
/**
		if (VATParameters::chimera)
		{
			// cout<<"Init chimeric alignment parameters"<<endl;
			// po::set_option(po::hit_cap, 15u);
			// po::set_option(po::min_identities, 14u);
			// po::set_option(po::padding, 8);
			// po::set_option(po::match, 1);
			// po::set_option(po::mismatch, -5);
			// po::set_option(po::lowmem, 1u);
			// cout<<VATParameters::hit_cap<<"\t"<<VATParameters::lowmem<<"\t"<<VATParameters::min_identities<<"\t"<<VATParameters::match<<"\t"<<VATParameters::mismatch<<endl;
			VATParameters::lowmem = 1u;
			VATParameters::penalty = -4;
			VATParameters::gap_extend = 0;
			VATParameters::gap_open = 0;
			VATParameters::match = 1;
			VATParameters::mismatch = -5;
			VATParameters::penalty = -4;
			VATParameters::padding = 8;
			VATParameters::hit_cap = 15u;
			VATParameters::min_identities = 14u;
		}
		else if (VATParameters::spilce || VATParameters::circ)
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
			VATParameters::gap_extend = 0;
			VATParameters::gap_open = 0;
			VATParameters::match = 5;
			VATParameters::mismatch = -2;
			// VATParameters::penalty = -4;
			VATParameters::padding = 8;
			VATParameters::hit_cap = 15;
			VATParameters::min_identities = 28;
		}
		else if (VATParameters::whole_genome_sequencing)
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
			VATParameters::lowmem = 1u;
			VATParameters::match = 5;
			VATParameters::mismatch = -2;
			VATParameters::penalty = -2;
			VATParameters::padding = 8;
			VATParameters::hit_cap = 2;
			VATParameters::min_identities = 28;
		}
		else if (VATParameters::whole_genome)
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
			VATParameters::padding = 8;
			VATParameters::hit_cap = 2;
			VATParameters::min_identities = 24;
			VATParameters::gapped_xdrop = 18;

		}
		else if (VATParameters::dna_homology)
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
