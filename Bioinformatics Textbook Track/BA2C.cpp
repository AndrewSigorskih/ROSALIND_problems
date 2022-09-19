#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>

using std::set;
using std::string;
using std::vector;
using matrix = vector<vector<double>>;
const char* INFILENAME = "rosalind_ba2c.txt";
const char* OUTFILENAME = "out.txt";

int getIndex(char const x)
{
    switch (x)
    {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return -1;
    }
}

class kmerProcessor
{
    int k;
    matrix profile;
    set<string> kmers;
public:
    kmerProcessor(std::ifstream& ist)
    {
        string text;
        getline(ist, text);
        ist >> k;
        for (int i = 0; i <= (text.size()-k); ++i)
            kmers.insert(text.substr(i, k));
        profile = matrix(4, vector<double>(k, 0.0));
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < k; ++j)
                ist >> profile[i][j];
    }
    string evaluateKmers();
};

string kmerProcessor::evaluateKmers()
{
    double max_prob = 0.0;
    string result;
    for (auto const& item: kmers)
    {
        double curprob = 1.0;
        for (int j = 0; j < k; ++j)
            curprob *= profile[getIndex(item[j])][j];
        if (curprob > max_prob)
        {
            max_prob = curprob;
            result = item;
        }
    }
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

    kmerProcessor processor(ist);
    ost << processor.evaluateKmers()<< "\n";
    return 0;
}