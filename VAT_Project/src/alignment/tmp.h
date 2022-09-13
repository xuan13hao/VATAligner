
void FindSeedChains(
    SFABuild &ref_seq, SFABuild &qry_seq, 
    AlignmentScoring &score, std::vector<SeedType> &seeds, 
    std::vector<SeedChainType> &seed_chains
) {

    if(seeds.size() <= 0)   return;
    // check the longest seed, if >90% of the query length, skip
    int qry_map_len = 0.9 * qry_seq.GetSeqLen(seeds[0].qry_id);
    for(int i = 0; i < seeds.size(); ++ i)   {
        if(seeds[i].qry_len >= qry_map_len) {
            SeedChainType sc;
            sc.chain.push_back(seeds[i]);
            sc.num_seeds = 1;
            seed_chains.push_back(sc);
        }
    }
    if(seed_chains.size() > 0) return;   
    
    // DEBUG
    //cout << "BEGIN of FindBestSplice:===" << endl;

    // sort the seeds based on locations
    sort(seeds.begin(), seeds.end(), SortSeedsByPos);
    // construct the pre-set for chaining, which contains the set of pre-seeds for each seed
    int **pre_seed_id = new int *[seeds.size()];   
    int **pre_gap_score = new int *[seeds.size()];      // assuming the seed break is a gap
    int **pre_spl_score = new int *[seeds.size()];      // assuming the seed break is a splice site
    for(int i = 0; i < seeds.size(); ++ i) {
        pre_seed_id[i] = new int [seeds.size()];
        pre_gap_score[i] = new int [seeds.size()];
        pre_spl_score[i] = new int [seeds.size()];
    }
    int *pre_count = new int [seeds.size()];    // stores the number of pre-seeds for each seed
    memset(pre_count, 0, seeds.size() * sizeof(int));
    
    // DEBUG
    //cout << "Finish initializing arrays" << endl;
    int seed_break, gap_size, overlap;
    score.GetMacroChain(seed_break, gap_size, overlap);
    Align align_obj;
    for(int i = 0; i < seeds.size() - 1; ++ i)   {
        int bound = qry_seq.GetSeqLen(seeds[0].qry_id);
        int j = i + 1;
        while(j < seeds.size()) {
            // set the right bound as the end position of the next seed
            if(seeds[i].qry_pos + seeds[i].qry_len - seeds[j].qry_pos <= overlap)    {
                bound = seeds[j].qry_pos + seeds[j].qry_len - 1;
                break;
            }
            ++ j;
        }
        // DEBUG
        //cout << "check index:\t" << i << "\t" << j << "\t" << bound << "\t" << seeds[j].qry_pos << "\t" << pre_count[j] << endl;
        
        // the remaining seeds must overlap with the seed for less than overlap, 
        // because all seeds are sorted by their starting locations
        /*
        int macro_gap_size_;    // the maximum allowed INDEL size (i.e. the difference in alignment length)
        int macro_seed_break_;  // the maximum seed break (i.e. the distance between two seeds)
        int macro_overlap_;     // the maximum seed overlap (i.e. the overlap between two seeds if they can be chained)*/
        for(; j < seeds.size(); ++ j)   {
            
            if(seeds[j].ref_id == seeds[i].ref_id && seeds[j].qry_pos <= bound)    {
                // here we need to consider splicing event and long-range gap
                int gap_score = -MAX;
                // if it fits the requirement of macro-gap
                int ref_gap = static_cast<int>(seeds[j].ref_pos) - static_cast<int>(seeds[i].ref_pos) - static_cast<int>(seeds[i].ref_len);
                int qry_gap = static_cast<int>(seeds[j].qry_pos) - static_cast<int>(seeds[i].qry_pos) - static_cast<int>(seeds[i].qry_len);
                // DEBUG
                //cout << "REF Seed POS:\t" << seeds[j].ref_pos << "\t" << seeds[i].ref_pos << "\t" << seeds[i].ref_len << endl;
                //cout << "QRY Seed POS:\t" << seeds[j].qry_pos << "\t" << seeds[i].qry_pos << "\t" << seeds[i].qry_len << endl;
                //cout << "GAP sizes:\t" << ref_gap << "\t" << qry_gap << endl;

                if(qry_gap <= seed_break && abs(ref_gap - qry_gap) <= gap_size)    
                {
                    int band = 2 * (abs(ref_gap - qry_gap) + 1);
                    gap_score = align_obj.AlignGlobalPairwise(
                        ref_seq.sequence_[seeds[i].ref_id] + seeds[i].ref_pos + 1, qry_seq.sequence_[seeds[i].qry_id] + seeds[i].qry_pos + 1,
                        ref_gap, qry_gap, band, score
                    );    
                }
                   
                int splice_score = IsSpliceJunction(ref_seq, qry_seq, seeds[i], seeds[j]) ? 0 : -MAX;
                // record the score
                if(gap_score > -MAX || splice_score > -MAX) {
                    // it means that the seed break might make sence, record
                    pre_seed_id[j][pre_count[j]] = i;
                    pre_gap_score[j][pre_count[j]] = gap_score;
                    pre_spl_score[j][pre_count[j]] = splice_score;
                    ++ pre_count[j];
                }

            }
        }
    }

    // use dynamic programming to find the best macro-gap or splice alignment
    int *p1_score = new int [seeds.size()]; int *p1_track = new int [seeds.size()]; bool *p1_splice = new bool [seeds.size()];
    // initialization
    p1_score[0] = seeds[0].score; p1_track[0] = -1;
    for(int i = 1; i < seeds.size(); ++ i)   {
        // perform dynamic programming, we need to record whether the seed break is a splicing site
        int max_score = seeds[i].score; int max_pre = i; bool is_splice = false;
        for(int j = 0; j < pre_count[i]; ++ j)   {
            int gs = pre_gap_score[i][j] + p1_score[pre_seed_id[i][j]];     // prev gap score
            int ss = pre_gap_score[i][j] + p1_score[pre_seed_id[i][j]];     // prev splice score
            if(gs > max_score)    {
                max_score = gs; max_pre = pre_seed_id[i][j]; is_splice = false;
            }
            if(ss > max_score)    {
                max_score = ss; max_pre = pre_seed_id[i][j]; is_splice = true;
            }
        }
        p1_score[i] = max_score;
        p1_track[i] = max_pre;
        p1_splice[i] = is_splice;
    }

    // traceback and recover the seed chains
    bool *tracked = new bool [seeds.size()];
    memset(tracked, false, seeds.size() * sizeof(bool));
    for(int j = seeds.size() - 1; j >= 0; j --)   {
        if(tracked[j]) continue;
        SeedChainType sc;
        int c = j;
        while(true) {
            tracked[c] = true;
            sc.chain.push_back(seeds[c]);
            sc.is_splice.push_back(p1_splice[c]);
            ++ sc.num_seeds;
            if(p1_track[c] == c)    {
                break;
            }   else
            {
                c = p1_track[c];
            }
        }
        // reverse the orders of the information in the seed chain
        std::reverse(sc.chain.begin(), sc.chain.end());
        std::reverse(sc.is_splice.begin(), sc.is_splice.end());
        seed_chains.push_back(sc);
    }

    // DEBUG
    //Seeding seed_obj;
    //for(int x = 0; x < seed_chains.size(); ++ x)   {
    //    cout << ">>>>>>>>>>>>>>>> Print sorted seeds: >>>>>>>>>>>>>>>>" << endl;
    //    for(int y = 0; y < seed_chains[x].chain.size(); ++ y)   {
    //        seed_obj.PrintSingleSeed(ref_seq, qry_seq, seed_chains[x].chain[y]);
    //    }
    //}

    // collect memory
    delete [] tracked;
    delete [] p1_score;
    delete [] p1_track;
    delete [] p1_splice;
    for(int i = 0; i < seeds.size(); ++ i)   {
        delete [] pre_seed_id[i];
        delete [] pre_gap_score[i];
        delete [] pre_spl_score[i];
    }
    delete [] pre_seed_id;
    delete [] pre_gap_score;
    delete [] pre_spl_score;
    delete [] pre_count;
    
    // DEBUG
    //cout << "DONE collectin memory" << endl;
    
    return;
}

/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2019 Jeremy Wang
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "klib/ksort.h"
#include "klib/khash.h"
#include "klib/kvec.h"
#include "chain.h"

#define score_pos_gt(a,b) ((a).score > (b).score)
KSORT_INIT(score_pos_cmp, score_pos, score_pos_gt)

#define pos_pair_lt(a,b) ((a).tpos < (b).tpos)
KSORT_INIT(pos_pair_cmp, posPair, pos_pair_lt)

chain* do_chain(khash_t(matchHash) *hits, int max_chains, int match_score, int max_gap, int min_chain_length) {

  // assess hits for each target
  pairVec anchors;
  khint_t bin;
  uint32_t target;

  scoreVec scores;
  kv_init(scores);        if (f[j] == )
        {
            /* code */
        }
        nto the scores vector of the current target - sets the backward limit for finding chained anchors

  // iterate through hits for each target, and append them to the same scores vector
  for (bin = kh_begin(hits); bin != kh_end(hits); ++bin) {
    if (!kh_exist(hits, bin)) continue;
    target = kh_key(hits, bin);
    anchors = kh_val(hits, bin);
    //fprintf(stderr, "Ref %u has %u anchors\n", target, kv_size(anchors));

    //sort anchor pairs by target pos increasing
    ks_mergesort(pos_pair_cmp, kv_size(anchors), anchors.a, 0);

    s.score = match_score;
    s.anchor_idx = 0;
    s.score_idx = ref_offset;
    s.prev = -1;
    s.used = 0;
    s.ref = target << 1 >> 1;
    s.rv = target >> 31;
    kv_push(score_pos, scores, s); // this actually copies the score_pos struct values, so we can reuse 's'

    int st = 0;

    for(i = 1; i < kv_size(anchors); i++) {

      //s.score = match_score;
      s.score = 1;
      s.anchor_idx = i;
      s.score_idx = i + ref_offset;
      s.prev = -1; // no predecessor

      /*
      while(st < i && kv_A(anchors, i).tpos - kv_A(anchors, st).tpos > 100000) st++;
      s.score = i - st;
      */

      for(j = i < h ? 0 : i-h; j < i; j++) {
        qdiff = (int)kv_A(anchors, i).qpos - (int)kv_A(anchors, j).qpos;
        tdiff = (int)kv_A(anchors, i).tpos - (int)kv_A(anchors, j).tpos;
        if(tdiff <= 0 || qdiff <= 0 || qdiff > max_gap || tdiff > max_gap) { // we often have multiple hits to the same target pos, so we can't let them chain with either other
          continue;
        }
        /*
        diffdiff = (qdiff > tdiff ? qdiff - tdiff : tdiff - qdiff);
        gap_cost = (diffdiff == 0 ? 0 : 0.01 * match_score * diffdiff + 0.5 * log2(diffdiff)); // affine gap penalty a la minimap2
        qdiff = (qdiff > tdiff ? tdiff : qdiff); // now not qdiff but the minimum difference
        score = kv_A(scores, j + ref_offset).score + (qdiff > match_score ? match_score : qdiff) - gap_cost;
        */
        score = kv_A(scores, j + ref_offset).score + 1; // scores are now just the # of properly ordered anchors
        if(score > s.score) {
          s.score = score;
          s.prev = j + ref_offset;
        }
      }

      kv_push(score_pos, scores, s);
    }
    ref_offset = kv_size(scores);
  }

  //fprintf(stderr, "%u scores\n", kv_size(scores));

  // copy unsorted scores:
  score_pos* anchor_scores = malloc(sizeof(score_pos) * kv_size(scores));
  memcpy(anchor_scores, scores.a, sizeof(score_pos) * kv_size(scores));
  //sort scores decreasing
  ks_mergesort(score_pos_cmp, kv_size(scores), scores.a, 0);

  // build non-overlapping chains from highest to lowest score
  chain* chains = calloc(max_chains, sizeof(chain));
  int c; // chain index
  i = 0; // scores index
  for(c = 0; c < max_chains && i < kv_size(scores); c++) {
    //fprintf(stderr, "building chain %d\n", c);
    int chain_len = 0;
    int chain_pos = kv_A(scores, i).score_idx; // get the pos which should be index into the anchor_scores array
    //fprintf(stderr, "chain pos: %d (of %u)\n", chain_pos, kv_size(scores));
    // get anchors associated with this score/ref
    bin = kh_get(matchHash, hits, kv_A(scores, i).ref + kv_A(scores, i).rv<<31);
    if(bin == kh_end(hits)) { // key not found, *shouldn't* happen
      fprintf(stderr, "something went very wrong with the chain computation!");
      return (chain*)NULL;
    }
    anchors = kh_val(hits, bin);
    //fprintf(stderr, "pos %d used: %u\n", chain_pos, anchor_scores[chain_pos].used);
    while(anchor_scores[chain_pos].used == 0) {
      //fprintf(stderr, "backtracking from %d (%d:%d) to %d\n", chain_pos, kv_A(anchors, anchor_scores[chain_pos].anchor_idx).qpos, kv_A(anchors, anchor_scores[chain_pos].anchor_idx).tpos, anchor_scores[chain_pos].prev);
      chain_len++;
      if(anchor_scores[chain_pos].prev == -1)
        break;
      chain_pos = anchor_scores[chain_pos].prev;
    }
    if(chain_len < min_chain_length) {
      c--;
      i++;
      continue;
    }
    kv_init(chains[c].anchors);
    chains[c].ref = kv_A(scores, i).ref;
    chains[c].rv = kv_A(scores, i).rv;
    chains[c].score = kv_A(scores, i).score;
    //fprintf(stderr, "chain %d (%c) (score %d, %d anchors) is from ref %u anchor %d\n", c, "+-"[kv_A(scores, i).rv], kv_A(scores, i).score, chain_len, chains[c].ref, kv_A(scores, i).anchor_idx);
    kv_resize(posPair, chains[c].anchors, chain_len);
    chains[c].anchors.n = chain_len; // set the size explicitly, then we'll set the values explicitly
    //fprintf(stderr, "creating chain %d of length %d\n", c, kv_size(chains[c].anchors));
    chain_pos = kv_A(scores, i).score_idx;
    for(j = 0; j < chain_len; j++) {
      anchor_scores[chain_pos].used = 1;
      kv_A(chains[c].anchors, chain_len-1-j) = kv_A(anchors, anchor_scores[chain_pos].anchor_idx);
      chain_pos = anchor_scores[chain_pos].prev;
    }
    i++;
  }
  //fprintf(stderr, "made %d chains\n", c);
  if(c < max_chains) kv_init(chains[c].anchors); // an empty vector will indicate the end of the chains array if there are fewer than max_chains results

  free(anchor_scores);
  kv_destroy(scores);
  return chains;
}

void free_chains(chain* chs) {
  int i;
  for(i = 0; ; i++) {
    if(kv_size(chs[i].anchors) == 0) {
      break;
    }
    kv_destroy(chs[i].anchors);
  }
  free(chs);
}