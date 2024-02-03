#ifndef __MINIMAP2_H__
#define __MINIMAP2_H__
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "../minimap2/minimap.h"
#include "../minimap2/kseq.h"
KSEQ_INIT(gzFile, gzread)

extern "C" {    

	struct VATMappedRec
	{
		std::string ref_name;
		int32_t ref_len;
		int32_t r_s;
		int32_t r_e;
		char strand;
		std::string qry_name;
		int32_t qry_len;
		int32_t q_s;
		int32_t q_e;
		int32_t score;
		std::string cigar_;
		int32_t mlen, blen;     // seeded exact match length; seeded alignment block length
		int32_t mapq;
	};
    struct CVector 
	{
        std::vector<VATMappedRec> cpp_vector;
    };

    CVector* create_vector() {
        return new CVector;
    }
    VATMappedRec get_element(CVector* cvec, int index) {
        // if (index < 0 || index >= cvec->cpp_vector.size()) {
        //     // Handle out-of-bounds error
        //     return nullptr; // You can choose a different error code
        // }
        return cvec->cpp_vector[index];
    }
	void push_back(CVector* cvec, VATMappedRec value) {
		cvec->cpp_vector.push_back(std::move(value));
	}

    int get_size(CVector* cvec) {
        return cvec->cpp_vector.size();
    }
    void destroy_vector(CVector* cvec) {
        delete cvec;
    }

	CVector* DNA_MM_module(const char* query_file, const char* subject_file)
	{
		mm_idxopt_t iopt;
		mm_mapopt_t mopt;
		// if ( VATParameters::short_model == 0 &&VATParameters::algn_type == VATParameters::dna)
		// {
		// 	mm_set_opt("short", &iopt, &mopt);
		// }else if(VATParameters::spilce  &&VATParameters::algn_type == VATParameters::dna )
		// {
		// 	mm_set_opt("splice", &iopt, &mopt);
		// }
		// // else if(VATParameters::spilce  &&VATParameters::algn_type == VATParameters::dna )
		// // {
		// // 	mm_set_opt("ava-ont", &iopt, &mopt);
		// // }
		// else
		// {
		mm_set_opt(0, &iopt, &mopt);
		// }
		int n_threads = 1;
		mm_verbose = 2; // disable message output to stderr
		// mm_set_opt(0, &iopt, &mopt);
		mopt.flag |= MM_F_CIGAR; // perform alignment
		// open query file for reading; you may use your favorite FASTA/Q parser
		gzFile f = gzopen(query_file, "r");
		assert(f);
		kseq_t *ks = kseq_init(f);
		// open index reader
		mm_idx_reader_t *r = mm_idx_reader_open(subject_file, &iopt, 0);
		mm_idx_t *mi;
		CVector* mr_v = create_vector();
		while ((mi = mm_idx_reader_read(r, n_threads)) != 0) { // traverse each part of the index
			mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
			mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
			gzrewind(f);
			kseq_rewind(ks);
			while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence & Read and print each sequence entry in the FASTA file
				mm_reg1_t *reg;
				int j, i, n_reg;
				reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0); // get all hits for the query
				for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
					VATMappedRec m_r;
					mm_reg1_t *r = &reg[j];
					assert(r->p); // with MM_F_CIGAR, this should not be NULL
					// printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
					m_r.ref_name = std::string(mi->seq[r->rid].name);
					m_r.r_s =r->rs;
					m_r.r_e =r->re;
					m_r.qry_len = ks->seq.l;
					m_r.strand = "+-"[r->rev];
					m_r.qry_name = std::string(ks->name.s);
					m_r.q_s = r->qs;
					m_r.q_e = r->qe;
					m_r.ref_len = mi->seq[r->rid].len;
					m_r.score = r->score;
					m_r.mlen = r->mlen;
					m_r.blen = r->blen;
					m_r.mapq = r->mapq;
					// printf("%s", m_r.ref_name);
					// printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
					// printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcg:Z:", mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re, r->mlen, r->blen, r->mapq);
					for (i = 0; i < r->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
					{
						// printf("%d%c", r->p->cigar[i]>>4, MM_CIGAR_STR[r->p->cigar[i]&0xf]);
						m_r.cigar_ = std::to_string(r->p->cigar[i]>>4)+MM_CIGAR_STR[r->p->cigar[i]&0xf];
					}
					push_back(mr_v,m_r);
					// putchar('\n');
					free(r->p);
				}
				free(reg);
			}
			mm_tbuf_destroy(tbuf);
			mm_idx_destroy(mi);
		}
		// int size = get_size(mr_v);
		// printf("Vector size: %d\n", size);
		// destroy_vector(mr_v);
		mm_idx_reader_close(r); // close the index reader
		kseq_destroy(ks); // close the query file
		gzclose(f);
		return mr_v;
		// return 0;
	}

void DNA_MultiThread_MM(const char* ref_file, const char** qry_file, int threads)
{

	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
    int ret;
	mm_set_opt(0, &iopt, &mopt);
	int threads = threads;
	mm_verbose = 0;
	mopt.flag |= MM_F_CIGAR; // perform alignment
	mm_idx_reader_t *r = mm_idx_reader_open(ref_file, &iopt, 0);
	mm_idx_t *mi;
	while ((mi = mm_idx_reader_read(r, threads)) != 0) { // traverse each part of the index
			mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
			ret = 0;
			//(const char**)&argv[2]
            ret = mm_map_file_frag(mi, 1, qry_file, &mopt, threads);
            if (ret < 0) break;
			mm_idx_destroy(mi);
		}
		mm_idx_reader_close(r); // close the index reader
}

const char** stringToConstCharPtrArray(const std::string& str) {
    // Allocate memory for a const char* array with one element
    const char** strArray = new const char*[1];

    // Create a copy of the string as a const char* and store it in the array
    strArray[0] = strdup(str.c_str());

    return strArray;
}



};
#endif // __MINIMAP2_H__