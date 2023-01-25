#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <set>
using namespace std;
using std::vector;
class DiagonalSeeds
{
    public:
	DiagonalSeeds() :
		len(0)
	{}
	DiagonalSeeds(int query_pos, int subject_pos, int len, int score) :
		i(query_pos),
		j(subject_pos),
		len(len),
		score(score)
	{
		
	}

	bool empty() const
	{
		return len == 0;
	}
	int subject_last() const
	{
		return j + len - 1;
	}
	int query_last() const
	{
		return i + len - 1;
	}
	int subject_end() const
	{
		return j + len;
	}
	int query_end() const
	{
		return i + len;
	}
	int diag() const
	{
		return i - j;
	}
	DiagonalSeeds transpose() const
	{
		return DiagonalSeeds(j, i, len, score);
	}
	int partial_score(int diff) const
	{
		return score*std::max(len - diff, 0) / len;
	}
	bool operator<=(const DiagonalSeeds &rhs) const
	{
		return i + len <= rhs.i && j + len <= rhs.j;
	}

	bool operator==(const DiagonalSeeds &rhs) const
	{
		return i == rhs.i && j == rhs.j && len == rhs.len;
	}
	static bool cmp_subject(const DiagonalSeeds &x, const DiagonalSeeds &y)
	{
		return x.j < y.j || (x.j == y.j && x.i < y.i);
	}
	static bool cmp_score(const DiagonalSeeds& x, const DiagonalSeeds& y)
	{
		return x.score > y.score;
	}
	static bool cmp_query_begin(const DiagonalSeeds &x, const DiagonalSeeds &y)
	{
		return x.i < y.i;
	}
	static bool cmp_subject_end(const DiagonalSeeds &x, const DiagonalSeeds &y)
	{
		return x.subject_end() < y.subject_end();
	}
	friend std::ostream& operator<<(std::ostream &s, const DiagonalSeeds &d)
	{
		s << "i=" << d.i << " j=" << d.j << " l=" << d.len << " score=" << d.score;
		return s;
	}
	int i, j, len, score;//query_pos, subject_pos
	int qry_id, sbj_id;
};


vector<DiagonalSeeds> chainingSeeds(vector<DiagonalSeeds>& diagonal_segment,int q_len,int max_gap)
{
    vector<DiagonalSeeds> chained_seed;
    int s_ =  diagonal_segment.size();
    int i, score, qdiff, tdiff, diffdiff, gap_cost, j, best_j;
    int h = 50; // number of previous anchors to check
    int match_score = 0;
    int *p1_score = new int [s_]; //chaining score
    int *p1_track = new int [s_]; // best predecessor

    p1_score[0] = diagonal_segment[0].score;
    p1_track[0] = -1;

    std::sort(diagonal_segment.begin(),diagonal_segment.end(),DiagonalSeeds::cmp_subject_end);
    for (size_t i = 1; i < diagonal_segment.size(); i++)
    {
        int tmp = diagonal_segment[i].score; 
        int max_pre = i;
        for(j = 0; j < i; j++) 
        {
            int tdiff = static_cast<int>(diagonal_segment[j].j) - static_cast<int>(diagonal_segment[i].j) - static_cast<int>(diagonal_segment[i].len);
            int qdiff = static_cast<int>(diagonal_segment[j].i) - static_cast<int>(diagonal_segment[i].i) - static_cast<int>(diagonal_segment[i].len);
            // cout<<"here"<<endl;
            // if(tdiff <= 0 || qdiff <= 0 || qdiff > max_gap || tdiff > max_gap) 
            // { // we often have multiple hits to the same target pos, so we can't let them chain with either other
            //     continue;
            // }
            // score = diagonal_segment[j].score;
            match_score = diagonal_segment[j].score;
            gap_cost = (diffdiff == 0 ? 0 : 0.01 * match_score * diffdiff + 0.5 * log2(diffdiff)); // affine gap penalty a la minimap2
            qdiff = (qdiff > tdiff ? tdiff : qdiff); // now not qdiff but the minimum difference
            score = diagonal_segment[j].score + (qdiff > match_score ? match_score : qdiff) - gap_cost; 
            // cout<<"score = "<<score<<endl;
             if (tmp < score)
             {  
                 tmp = score;
                 max_pre = j;
             } 
        }

        p1_score[i] = tmp;
        p1_track[i] = max_pre;

    }

    for (size_t i = 0; i < s_; i++)
    {
        cout<<p1_track[i]<<"\t"<<p1_score[i]<<endl;
    }
    int best;
    set<int> visited;
    set<int>::iterator it;
    for (size_t i = 0; i < s_; i++)
    {
        visited.insert(p1_track[i]);
    }
    for (it = visited.begin(); it != visited.end(); ++it)
    {
        cout<<*it<<endl;
        if (*it == -1)
        {
            chained_seed.push_back(diagonal_segment[0]);
        }else
        {
            chained_seed.push_back(diagonal_segment[*it]);
        }   
        
    }
    cout<<chained_seed.size()<<endl;
    
    /*
    for (size_t j = s_ - 1; j >= 0; --j)
    { 
        if (visited.count(j) > 0)
        {
            continue;
        }
        int current = j;
        
        while (p1_track[current] != -1 || visited.count(p1_track[current]) < 0)
        {
            cout<<current<<"\t"<<p1_track[current]<<endl;
            visited.insert(current);
            chained_seed.push_back(diagonal_segment[current]);
            current = p1_track[current];
        
        }
    }
    */
    // std::reverse(chained_seed.begin(),chained_seed.end());
    delete [] p1_score;
    delete [] p1_track;

    return chained_seed;
}



int main()
{
    DiagonalSeeds ds1(1,2,10,15);
    DiagonalSeeds ds2(1,2,10,40);//*
    DiagonalSeeds ds3(1,2,10,15);//*
    DiagonalSeeds ds4(100,103,13,40);
    DiagonalSeeds ds5(1,1,15,15);//*
    // DiagonalSeeds ds6(49,52,8,8);//*
    vector<DiagonalSeeds> vds;
    vector<DiagonalSeeds> result;
    // vector<SeedChainType> seed_chains;

    vds.push_back(ds1);
    vds.push_back(ds2);
    vds.push_back(ds3);
    vds.push_back(ds4);
    vds.push_back(ds5);
    // vds.push_back(ds6);
    // cout<<"size = "<<vds.size()<<endl;
    result = chainingSeeds(vds,50,10);
    // result = findOptimalSeeds(vds,30,15);
    // findSeedChain(vds,seed_chains,30);
    // cout<<
    cout<<"result = "<<result.size()<<endl;

    for (size_t i = 0; i < result.size(); i++)
    {
        cout<<result[i].i<<"\t"<<result[i].j<<"\t"<<result[i].len<<"\t"<<result[i].score<<endl;
    }
    
    return 1;
}

