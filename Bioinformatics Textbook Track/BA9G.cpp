#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm> // sort

using std::vector;
using std::string;
const char* INFILENAME = "rosalind_ba9g.txt";
const char* OUTFILENAME = "out.txt";

// time complexity O(n*logn*logn) implementation
// source: https://www.geeksforgeeks.org/suffix-array-set-2-a-nlognlogn-algorithm/

struct suffix
{
    int index;
    int rank[2];
    bool operator < (const suffix& other)
    {
        return (this->rank[0] == other.rank[0]) ? (this->rank[1] < other.rank[1] ? 1: 0) :
               (this->rank[0] < other.rank[0] ? 1: 0);
    }
};

vector<int> buildSuffixArray(string const& text)
{
    int n = text.size();
    suffix suffixes[n];

    for (int i = 0; i < n; ++i)
    {
        suffixes[i].index = i;
        suffixes[i].rank[0] = text[i] - 'a';
        suffixes[i].rank[1] = ((i+1) < n) ? (text[i+1] - 'a') : -1;
    }

    std::sort(suffixes, suffixes+n);

    int ind[n];
    for (int k = 4; k < 2*n; k = k*2)
    {
        int rank = 0;
        int prev_rank = suffixes[0].rank[0];
        suffixes[0].rank[0] = rank;
        ind[suffixes[0].index] = 0;

        for (int i = 1; i < n; ++i)
        {
            if (suffixes[i].rank[0] == prev_rank &&
                    suffixes[i].rank[1] == suffixes[i-1].rank[1])
            {
                prev_rank = suffixes[i].rank[0];
                suffixes[i].rank[0] = rank;
            } else {
                prev_rank = suffixes[i].rank[0];
                suffixes[i].rank[0] = ++rank;
            }
            ind[suffixes[i].index] = i;
        }

        for (int i = 0; i < n; i++)
        {
            int nextindex = suffixes[i].index + k/2;
            suffixes[i].rank[1] = (nextindex < n) ?
                                  suffixes[ind[nextindex]].rank[0] : -1;
        }
        std::sort(suffixes, suffixes+n);
    }

    vector<int> result(n);
    for (int i = 0; i < n; i++)
        result[i] = suffixes[i].index;
    return result;
}


int main()
{
    std::ifstream ist{INFILENAME};
    if (!ist) 
    {
        std::cerr << "Cannot open input file!\n";
        exit(1);
    }
    std::ofstream ost{OUTFILENAME};
    if (!ost) 
    {
        std::cerr << "Cannot open output file!\n";
        exit(1);
    }

    string input;
    getline(ist, input);

    auto result = buildSuffixArray(input);
    for (int i = 0; i < result.size() - 1; ++i)
        ost << result[i] << ", ";
    ost << result[result.size()-1] << "\n";
}