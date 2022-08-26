#include <iostream>
#include "../basic/DiagonalSeeds.h"
#include "./FindSeedsChain.h"
using namespace std;
//vector<DiagonalSeeds> findOptimalSeeds(vector<DiagonalSeeds>& diagonal_segment,int q_len,int max_gap)


int main()
{
    DiagonalSeeds ds1(1,2,10,5);
    DiagonalSeeds ds2(7,9,9,5);//*
    DiagonalSeeds ds3(10,12,11,5);//*
    DiagonalSeeds ds4(17,19,13,5);
    DiagonalSeeds ds5(21,22,15,5);//*
    DiagonalSeeds ds6(27,29,8,5);//*
    DiagonalSeeds ds7(31,42,12,5);//*
    DiagonalSeeds ds8(37,49,9,5);//*
    vector<DiagonalSeeds> vds;
    vector<DiagonalSeeds> result;
    vds.push_back(ds1);
    vds.push_back(ds2);
    vds.push_back(ds3);
    vds.push_back(ds4);
    vds.push_back(ds5);
    vds.push_back(ds6);
    vds.push_back(ds7);
    vds.push_back(ds8);
    cout<<"size = "<<vds.size()<<endl;
    result = findOptimalSeeds(vds,30,10);
    cout<<"result = "<<result.size()<<endl;

    for (size_t i = 0; i < result.size(); i++)
    {
        cout<<"len = "<<result[i].len<<endl;
    }
    
    return 1;
}