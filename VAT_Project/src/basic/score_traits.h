
#ifndef SCORE_TRAITS_H_
#define SCORE_TRAITS_H_

template<typename _val>
uint8_t blast_seq_code()
{ return 0; }

template<>
uint8_t blast_seq_code<Protein>()
{ return BLASTAA_SEQ_CODE; }

template<typename _val>
int16_t blast_load_karlin_blk(Blast_KarlinBlk* kbp,
		Blast_KarlinBlk* kbp_ungap,
		int gap_open,
		int gap_extend,
		int reward,
		int penalty,
		const char *matrix)
{ return 0; }

template<>
int16_t blast_load_karlin_blk<Protein>(Blast_KarlinBlk* kbp,
		Blast_KarlinBlk* kbp_ungap,
		int gap_open,
		int gap_extend,
		int reward,
		int penalty,
		const char *matrix)
{
	return Blast_KarlinBlkGappedLoadFromTables(kbp,
			gap_open,
			gap_extend,
			matrix);
}

template<typename _val>
const uint8_t* blast_alphabet()
{ return 0; }

template<>
const uint8_t* blast_alphabet<Protein>()
{ return AMINOACID_TO_NCBISTDAA; }

#ifdef EXTRA
#include "../../../extra/score_traits.h"
#endif

#endif /* SCORE_TRAITS_H_ */