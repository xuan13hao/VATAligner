#ifndef __REPORT_H__
#define __REPORT_H__
#include <iostream>
#include <fstream>
#include "../alignment/FindSeedsChain.h"
using namespace std;
bool SortAlnByScore(const ChimeraAlnType &a, const ChimeraAlnType &b)  {
  if(a.score > b.score || (a.score == b.score && a.arm1.sbj_id > b.arm1.sbj_id))  {
    return true;
  }  
  return false;
}
void ReportBLASTFormat(
    std::vector<ChimeraAlnType> &aln, const int &num_best, 
    std::ofstream &out_fh ) {
    if(aln.size() <= 0) return;

    // sort the alignment segments by score
    sort(aln.begin(), aln.end(), SortAlnByScore);
    
    // TODO:
    // this is a simplified version, we need to get the actual name ID in the future
    int prev_score = aln[0].score;   
    int num_to_print = num_best;
    for(int i = 0; i < aln.size(); ++ i)   {
        if(aln[i].score < prev_score)   {
            prev_score = aln[i].score;
            -- num_to_print;
        }
        if(num_to_print <= 0)   break;
        // print the first arm (no matter whether the alignment is chimera)
        out_fh << aln[i].arm1.qry_id << "\t" << aln[i].arm1.sbj_id << "\t";
        out_fh << aln[i].arm1.i << "\t" << aln[i].arm1.i + aln[i].arm1.len - 1 << "\t";
        out_fh << aln[i].arm1.j << "\t" << aln[i].arm1.j + aln[i].arm1.len - 1 << "\t";
        out_fh << aln[i].score << "\t" << i << endl;
        if(aln[i].is_chimera)    {
            // if it is chimera, print the second
            out_fh << aln[i].arm2.qry_id << "\t" << aln[i].arm2.sbj_id << "\t";
            out_fh << aln[i].arm2.i << "\t" << aln[i].arm2.i + aln[i].arm2.len - 1 << "\t";
            out_fh << aln[i].arm2.j << "\t" << aln[i].arm2.j + aln[i].arm2.len - 1 << "\t";
            out_fh << aln[i].score << "\t" << i << endl;
        }  
    }
    return;
}

#endif // __REPORT_H__