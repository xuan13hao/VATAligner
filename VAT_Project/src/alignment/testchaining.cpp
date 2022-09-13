#include <iostream>
#include "../basic/DiagonalSeeds.h"
#include "./FindSeedsChain.h"
using namespace std;
//vector<DiagonalSeeds> findOptimalSeeds(vector<DiagonalSeeds>& diagonal_segment,int q_len,int max_gap)
//void findSeedChain(vector<DiagonalSeeds>& diagonal_segment,
//vector<SeedChainType>& seed_chains,int q_len)


int main()
{
    DiagonalSeeds ds1(1,2,10,15);
    DiagonalSeeds ds2(11,12,9,16);//*
    DiagonalSeeds ds3(21,22,11,18);//*
    DiagonalSeeds ds4(31,32,13,15);
    DiagonalSeeds ds5(41,42,15,14);//*
    DiagonalSeeds ds6(51,52,8,16);//*
    vector<DiagonalSeeds> vds;
    vector<DiagonalSeeds> result;
    vector<SeedChainType> seed_chains;

    vds.push_back(ds1);
    vds.push_back(ds2);
    vds.push_back(ds3);
    vds.push_back(ds4);
    vds.push_back(ds5);
    vds.push_back(ds6);
    cout<<"size = "<<vds.size()<<endl;
    result = findOptimalSeeds(vds,30,15);
    // findSeedChain(vds,seed_chains,30);
    // cout<<
    cout<<"result = "<<result.size()<<endl;

    for (size_t i = 0; i < result.size(); i++)
    {
        cout<<result[i].i<<"\t"<<result[i].j<<"\t"<<result[i].len<<"\t"<<result[i].score<<endl;
    }
    
    return 1;
}