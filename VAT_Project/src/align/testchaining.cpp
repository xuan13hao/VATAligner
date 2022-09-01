#include <iostream>
#include "../basic/DiagonalSeeds.h"
#include "./FindSeedsChain.h"
using namespace std;
//vector<DiagonalSeeds> findOptimalSeeds(vector<DiagonalSeeds>& diagonal_segment,int q_len,int max_gap)
//void findSeedChain(vector<DiagonalSeeds>& diagonal_segment,
//vector<SeedChainType>& seed_chains,int q_len)


int main()
{
    DiagonalSeeds ds1(1,5,10,15);
    DiagonalSeeds ds2(17,29,9,16);//*
    DiagonalSeeds ds3(1,32,11,18);//*
    DiagonalSeeds ds4(17,19,13,15);
    DiagonalSeeds ds5(11,22,15,14);//*
    DiagonalSeeds ds6(17,29,8,16);//*
    DiagonalSeeds ds7(31,42,12,10);//*
    DiagonalSeeds ds8(37,49,9,12);//*
    vector<DiagonalSeeds> vds;
    vector<DiagonalSeeds> result;
    vector<SeedChainType> seed_chains;

    vds.push_back(ds1);
    vds.push_back(ds2);
    vds.push_back(ds3);
    vds.push_back(ds4);
    vds.push_back(ds5);
    vds.push_back(ds6);
    vds.push_back(ds7);
    vds.push_back(ds8);
    cout<<"size = "<<vds.size()<<endl;
    result = findOptimalSeeds(vds,30,5);
    // findSeedChain(vds,seed_chains,30);
    // cout<<
    cout<<"result = "<<result.size()<<endl;

    for (size_t i = 0; i < result.size(); i++)
    {
        cout<<"len = "<<result[i].len<<endl;
    }
    
    return 1;
}