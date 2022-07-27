#include <seqan/align.h>
#include <seqan/stream.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>
#include <iostream>
#include "basic/Minimizer.h"
using namespace std;
using namespace seqan;

int main()
{
    // typedef Seed<Simple>    TSeed;
    // typedef SeedSet<TSeed> TSeedSet;

    // TSeedSet seedSet;
    // addSeed(seedSet, TSeed(0, 0, 2), Single());
    // addSeed(seedSet, TSeed(3, 5, 2), Single());
    // addSeed(seedSet, TSeed(4, 2, 3), Single());
    // addSeed(seedSet, TSeed(9, 9, 2), Single());

    // std::cout << "Resulting seeds.\n";
    // typedef Iterator<TSeedSet>::Type TIter;
    // for (TIter it = begin(seedSet, Standard()); it != end(seedSet, Standard()); ++it)
    //     std::cout << "(" << beginPositionH(*it) << ", " << endPositionH(*it)
    //               << ", " << beginPositionV(*it) << ", " << endPositionV(*it)
    //               << ", " << lowerDiagonal(*it) << ", " << upperDiagonal(*it)
    //               << ")\n";
/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param p      minimizers; p->a[i].x is the 2k-bit hash value;
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
    const char *str = "GTAAACCTTAAT";
    mm128_v *p;
    p = (mm128_v*)calloc(1, sizeof(mm128_v));
    mm_sketch(str, 4, 4, 3, 1, p);
    for (size_t i = 0; i < p->n; i++)
    {
        // cout<<"x = "<<p->m<<",y = "<<p->n<<endl;
        cout<<"x = "<<p->a[i].x<<",y = "<<p->a[i].y<<endl;
    }
    
    
    return 0;
}