#ifndef __PAIRCHIMERA_H__
#define __PAIRCHIMERA_H__
#include <vector>
#include <algorithm> 
#include <string.h>
#include <math.h>
#include <assert.h>
#include "../basic/Seed.h"
template<typename _locr, typename _locl>
struct ChimeraAlnType  {
    DiagonalSeeds<_locr, _locl> arm1, arm2;
    bool is_chimera;
    int score;

    ChimeraAlnType& operator =(const ChimeraAlnType& a)
    {
        arm1 = a.arm1;
        arm2 = a.arm2;
        is_chimera = a.is_chimera;
        score = a.score;
        return *this;
    }
};

template<typename _locr, typename _locl>
void pairChimeraSeeds(std::vector<DiagonalSeeds<_locr, _locl> > &seeds, 
std::vector<ChimeraAlnType<_locr, _locl> > &paired_seeds,int q_len) {
    if(seeds.size() == 0)   return;
    if(seeds.size() == 1)   
    {
        // if it contains only one seed, take it directly
        ChimeraAlnType<_locr, _locl> p;
        p.arm1 = seeds[0]; p.score = seeds[0].score; p.is_chimera = false;
        paired_seeds.push_back(p);
        return;
    }
    // check whether all seeds are from the same query
    int q_id = seeds[0].qry_id;
    for(int i = 1; i < seeds.size(); ++ i)   
    {
        assert(q_id == seeds[i].qry_id);
    }
    // sort the seeds based on their query positions
    std::sort(seeds.begin(),seeds.end(),DiagonalSeeds<_locr, _locl>::cmp_query_begin);
    int q_seq_len = q_len;
    // int q_seq_len = qry_seq.GetSeqLen(q_id);
    int **cand = new int* [q_seq_len];
    for(int i = 0; i < q_seq_len; ++ i)   
    {
        cand[i] = new int [seeds.size()];
    }
    int *cand_count = new int [q_seq_len];
    memset(cand_count, 0, q_seq_len * sizeof(int));
    for(int j = 0; j < seeds.size(); ++ j)   
    {
        int p = seeds[j].i;
        cand[p][cand_count[p] ++] = j;
    }
    // for(int x = 0; x < q_seq_len; ++ x)   
    // {
    //    cout << "Num candidates:    " << cand_count[x] << endl;
    //    for(int y = 0; y < cand_count[x]; ++ y)   {
    //        cout << cand[x][y] << endl;
    //    }
    // }

    int c_pen = -5; 
    int c_overlap = 5;
    bool *taken = new bool [seeds.size()];
    memset(taken, false, seeds.size() * sizeof(bool));
    // find all chimera alignment that improves the singular mapping score
    for(int i = 0; i < seeds.size(); ++ i)   
    {
        int bound = seeds[i].i + seeds[i].len - c_overlap;
        assert(bound < q_seq_len);
        int max_score = seeds[i].score;
        ChimeraAlnType<_locr, _locl> chim;
        chim.arm1 = seeds[i]; chim.is_chimera = false; chim.score = seeds[i].score;
        // cout<<chim.arm1.qry_id<<"\t"<<chim.arm1.sbj_id<<"\t"<<chim.arm1.i<<"\t"<<chim.arm1.j<<"\t"<<chim.score <<endl;
        // cout<<seeds[i].qry_id<<"\t"<<seeds[i].sbj_id<<"\t"<<seeds[i].i<<"\t"<<seeds[i].j<<endl;
        vector<int> p_seeds;
        for(int j = bound; j < q_seq_len; ++ j)   
        {
            for(int k = 0; k < cand_count[j]; ++ k)   
            {
                if(seeds[cand[j][k]].score + c_pen > 0)    
                {
                    if(seeds[i].score + seeds[cand[j][k]].score + c_pen >= max_score)    
                    {
                        
                        max_score = seeds[i].score + seeds[cand[j][k]].score + c_pen;
                        // cout<<"cand = "<<cand[j][k]<<endl;
                        p_seeds.push_back(cand[j][k]);
                    }
                }
            }
        }
        // cout<<"p_seeds.size() = "<<p_seeds.size()<<endl;
        if(p_seeds.size() > 0)   
        {
            for(int l = 0; l < p_seeds.size(); ++ l)   
            {
                if(seeds[i].score + seeds[p_seeds[l]].score + c_pen >= max_score)    
                {
                    // a chimeric mapping is detected
                    ChimeraAlnType<_locr, _locl> cm = chim;
                    cm.arm2 = seeds[p_seeds[l]];
                    cm.is_chimera = true; cm.score = max_score;
                    // cout<<chim.score<<"\t"<<chim.is_chimera<<"\t"<<chim.arm1.i<<"\t"<<chim.arm1.j<<endl;

                    taken[i] = taken[p_seeds[l]] = true;
                    paired_seeds.push_back(cm);
                }
            }
        }  
        else    
        {
            // no chimeric mapping can be detected
            if(!taken[i])    
            {
                // only include it in the final set if the alignment has not been taken 
                // cout<<chim.arm1.qry_id<<"\t"<<chim.arm1.sbj_id<<"\t"<<chim.arm1.i<<"\t"<<chim.arm1.j<<endl;
                paired_seeds.push_back(chim);
            }
            taken[i] = true;
        }
    }

    // collect memory
    for(int i = 0; i < q_seq_len; ++ i)   
    {
        delete [] cand[i];
    }
    delete [] cand;
    delete [] cand_count;

    return;
}

#endif // __PAIRCHIMERA_H__