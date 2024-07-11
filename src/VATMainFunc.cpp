#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
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
using std::string;
using std::vector;
using std::map;

class ArgumentParser {
public:
    ArgumentParser(int &argc, char **argv) {
        for (int i = 1; i < argc; ++i) {
            tokens.push_back(string(argv[i]));
        }
    }

    const string& getCmdOption(const string &option) const {
        vector<string>::const_iterator itr;
        itr = std::find(tokens.begin(), tokens.end(), option);
        if (itr != tokens.end() && ++itr != tokens.end()) {
            return *itr;
        }
        static const string empty_string("");
        return empty_string;
    }

    bool cmdOptionExists(const string &option) const {
        return std::find(tokens.begin(), tokens.end(), option) != tokens.end();
    }

private:
    vector<string> tokens;
};

void setDefaultValues() {
    VATParameters::threads_ = 1;
    VATParameters::max_alignments = 25;
    VATParameters::toppercent = 98.0;
    VATParameters::max_evalue = 0.001;
    VATParameters::min_bit_score = 0.0;
    VATParameters::min_id = 0.0;
    VATParameters::tmpdir = "/dev/shm";
    VATParameters::gap_open = -1;
    VATParameters::gap_extend = -1;
    VATParameters::reward = 2;
    VATParameters::seed_len = 15;
    VATParameters::penalty = -3;
    VATParameters::match = 5;
    VATParameters::mismatch = -4;
    VATParameters::spaced_seed = "null";
    VATParameters::matrix = "blosum62";
    VATParameters::max_seed_freq = -20.0;
    VATParameters::hit_cap = 0;
    VATParameters::min_identities = 0;
    VATParameters::xdrop = 18;
    VATParameters::gapped_xdrop = 18;
    VATParameters::min_ungapped_raw_score = 0;
    VATParameters::min_hit_score = 0;
    VATParameters::padding = 0;
    VATParameters::output_format = "tab";
}

int main(int argc, char* argv[]) {
    ArgumentParser input(argc, argv);
    string command;

    try {
        setDefaultValues();

        if (input.cmdOptionExists("-h") || input.cmdOptionExists("--help")) {
            cout << endl << "Syntax:" << endl;
            cout << "  VAT COMMAND [OPTIONS]" << endl << endl;
            cout << "Commands:" << endl;
            cout << "makevatdb\tBuild VAT database from a FASTA file" << endl;
            cout << "protein\tAlign protein query sequences against a protein reference database" << endl;
            cout << "dna\tAlign DNA query sequences against a DNA reference database" << endl;
            cout << "view\tView VAT alignment archive (vaa) formatted file" << endl;
            cout << endl;
            cout << "General options:" << endl;
            cout << "  --help, -h        produce help message" << endl;
            cout << "  --threads, -p     number of cpu threads (default=1)" << endl;
            cout << "  --db, -d          database file" << endl;
            cout << "  --vaa, -a         VAT alignment archive (vatr) file" << endl;
            cout << "  --dbtype          database type (nucl/prot)" << endl;
            cout << "Makedb options:" << endl;
            cout << "  --in, -i          input reference file in FASTA format" << endl;
            cout << "Aligner options:" << endl;
            cout << "  --query, -q       input query file" << endl;
            cout << "  --maxtarget_seqs, -k   maximum number of target sequences to report alignments for (default=25)" << endl;
            cout << "  --top             report alignments within this percentage range of top alignment score (default=98)" << endl;
            cout << "  --evalue, -e      maximum e-value to report alignments (default=0.001)" << endl;
            cout << "  --min_score       minimum bit score to report alignments (default=0)" << endl;
            cout << "  --report_id       minimum identity% to report an alignment (default=0)" << endl;
            cout << "  --tmpdir, -t      directory for temporary files (default=/dev/shm)" << endl;
            cout << "  --gapopen         gap open penalty, -1=default (11 for protein)" << endl;
            cout << "  --gapextend       gap extension penalty, -1=default (1 for protein)" << endl;
            cout << "  --reward          match reward score (blastn only, default=2)" << endl;
            cout << "  --seed_len, -S    Seed length (default=15 for DNA, 8 for Protein)" << endl;
            cout << "  --penalty         mismatch penalty score (blastn only, default=-3)" << endl;
            cout << "  --match           match score (default=5)" << endl;
            cout << "  --mismatch        mismatch score (default=-4)" << endl;
            cout << "  --simd_sort       Double-index based on SIMD" << endl;
            cout << "  --chimera         Chimera alignment" << endl;
            cout << "  --circ            Circ alignment" << endl;
            cout << "  --wga             Whole-genome alignment" << endl;
            cout << "  --wgs             Whole-genome sequencing" << endl;
            cout << "  --splice          Splice alignments" << endl;
            cout << "  --dnah            DNA homology" << endl;
            cout << "  --spaced          Spaced seed (default=null)" << endl;
            cout << "  --matrix          score matrix for protein alignment (default=blosum62)" << endl;
            cout << "Advanced options:" << endl;
            cout << "  --seed_freq       maximum seed frequency (default=-20)" << endl;
            cout << "  --max_hits, -C    maximum number of hits to consider for one seed (default=0)" << endl;
            cout << "  --pre_filter      minimum number of identities for pre-filter hit (default=0)" << endl;
            cout << "  --xdrop           xdrop for ungapped alignment (default=18)" << endl;
            cout << "  --gapped_xdrop, -X  xdrop for gapped alignment in bits (default=18)" << endl;
            cout << "  --ungapped_score  minimum raw alignment score to continue local extension (default=0)" << endl;
            cout << "  --hit_score       minimum score to keep a tentative alignment (default=0)" << endl;
            cout << "  --band            band for dynamic programming computation (default=0)" << endl;
            cout << "  --for_only        only forward strand alignment" << endl;
            cout << "View options:" << endl;
            cout << "  --out, -o         output file" << endl;
            cout << "  --outfmt, -f      output format (default=tab)" << endl;
            return 0;
        }

        if (input.cmdOptionExists("--threads")) VATParameters::threads_ = std::stoi(input.getCmdOption("--threads"));
        if (input.cmdOptionExists("--db")) VATParameters::database = input.getCmdOption("--db");
        if (input.cmdOptionExists("--vaa")) VATParameters::daa_file = input.getCmdOption("--vaa");
        if (input.cmdOptionExists("--dbtype")) VATParameters::db_type = input.getCmdOption("--dbtype");
        if (input.cmdOptionExists("--in")) VATParameters::input_ref_file = input.getCmdOption("--in");
        if (input.cmdOptionExists("--query")) VATParameters::query_file = input.getCmdOption("--query");
        if (input.cmdOptionExists("--maxtarget_seqs")) VATParameters::max_alignments = std::stoull(input.getCmdOption("--maxtarget_seqs"));
        if (input.cmdOptionExists("--top")) VATParameters::toppercent = std::stod(input.getCmdOption("--top"));
        if (input.cmdOptionExists("--evalue")) VATParameters::max_evalue = std::stod(input.getCmdOption("--evalue"));
        if (input.cmdOptionExists("--min_score")) VATParameters::min_bit_score = std::stod(input.getCmdOption("--min_score"));
        if (input.cmdOptionExists("--report_id")) VATParameters::min_id = std::stod(input.getCmdOption("--report_id"));
        if (input.cmdOptionExists("--tmpdir")) VATParameters::tmpdir = input.getCmdOption("--tmpdir");
        if (input.cmdOptionExists("--gapopen")) VATParameters::gap_open = std::stoi(input.getCmdOption("--gapopen"));
        if (input.cmdOptionExists("--gapextend")) VATParameters::gap_extend = std::stoi(input.getCmdOption("--gapextend"));
        if (input.cmdOptionExists("--reward")) VATParameters::reward = std::stoi(input.getCmdOption("--reward"));
        if (input.cmdOptionExists("--seed_len")) VATParameters::seed_len = std::stoi(input.getCmdOption("--seed_len"));
        if (input.cmdOptionExists("--penalty")) VATParameters::penalty = std::stoi(input.getCmdOption("--penalty"));
        if (input.cmdOptionExists("--match")) VATParameters::match = std::stoi(input.getCmdOption("--match"));
        if (input.cmdOptionExists("--mismatch")) VATParameters::mismatch = std::stoi(input.getCmdOption("--mismatch"));
        if (input.cmdOptionExists("--spaced")) VATParameters::spaced_seed = input.getCmdOption("--spaced");
        if (input.cmdOptionExists("--matrix")) VATParameters::matrix = input.getCmdOption("--matrix");
        if (input.cmdOptionExists("--seed_freq")) VATParameters::max_seed_freq = std::stod(input.getCmdOption("--seed_freq"));
        if (input.cmdOptionExists("--max_hits")) VATParameters::hit_cap = std::stoul(input.getCmdOption("--max_hits"));
        if (input.cmdOptionExists("--pre_filter")) VATParameters::min_identities = std::stoul(input.getCmdOption("--pre_filter"));
        if (input.cmdOptionExists("--xdrop")) VATParameters::xdrop = std::stoi(input.getCmdOption("--xdrop"));
        if (input.cmdOptionExists("--gapped_xdrop")) VATParameters::gapped_xdrop = std::stoi(input.getCmdOption("--gapped_xdrop"));
        if (input.cmdOptionExists("--ungapped_score")) VATParameters::min_ungapped_raw_score = std::stoi(input.getCmdOption("--ungapped_score"));
        if (input.cmdOptionExists("--hit_score")) VATParameters::min_hit_score = std::stoi(input.getCmdOption("--hit_score"));
        if (input.cmdOptionExists("--band")) VATParameters::padding = std::stoi(input.getCmdOption("--band"));
        if (input.cmdOptionExists("--out")) VATParameters::output_file = input.getCmdOption("--out");
        if (input.cmdOptionExists("--outfmt")) VATParameters::output_format = input.getCmdOption("--outfmt");

        VATParameters::forward_only = input.cmdOptionExists("--for_only");
        VATParameters::chimera = input.cmdOptionExists("--chimera");
        VATParameters::whole_genome_sequencing = input.cmdOptionExists("--wgs");
        VATParameters::dna_homology = input.cmdOptionExists("--dnah");
        VATParameters::whole_genome = input.cmdOptionExists("--wga");
        VATParameters::circ = input.cmdOptionExists("--circ");
        VATParameters::spilce = input.cmdOptionExists("--splice");

        VATParameters::alignment_traceback = !input.cmdOptionExists("--no-traceback");

        if (input.cmdOptionExists("makevatdb") && !VATParameters::input_ref_file.empty() && !VATParameters::database.empty()) {
            if (VATParameters::db_type == "nucl") {
                VATParameters::chunk_size = 8;
                RunModel::CreateDNADB();
            } else if (VATParameters::db_type == "prot") {
                VATParameters::chunk_size = 8;
                RunModel::CreateProteinDB();
            } else {
                cerr << "Failed to get database type. Please refer to the readme for instructions." << endl;
            }
        } else if ((input.cmdOptionExists("protein") || input.cmdOptionExists("dna") || input.cmdOptionExists("blastx")) &&
                   !VATParameters::query_file.empty() && !VATParameters::database.empty() && !VATParameters::daa_file.empty()) {
            if (input.cmdOptionExists("protein")) {
                RunModel::ProteinAlign();
            } else if (input.cmdOptionExists("dna")) {
                RunModel::DNAAlign();
            } else if (input.cmdOptionExists("blastx")) {
                RunModel::BlastxAlign();
            } else {
                cerr << "Failed to get alignment type. Please refer to the readme for instructions." << endl;
            }
        } else if (input.cmdOptionExists("view") && !VATParameters::daa_file.empty()) {
            view();
        } else {
            cout << "Insufficient arguments. Use VAT -h for help.\n";
        }
    } catch (std::bad_alloc &e) {
        cerr << "Failed to allocate sufficient memory. Please refer to the readme for instructions on memory usage." << endl;
    } catch (std::exception &e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    } catch (...) {
        cerr << "Exception of unknown type!\n";
        return 1;
    }

    return 0;
}
