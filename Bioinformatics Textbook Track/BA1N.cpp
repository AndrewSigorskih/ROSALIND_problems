#include <iostream>
#include <fstream>
#include <string>
// naive approach: just check all possible kmers
using std::string;
const char* INFILENAME = "rosalind_ba1n.txt";
const char* OUTFILENAME = "out.txt";

bool hammingDistanceLEd(string const& s1, string const& s2, int const d)
{
    if (s1.size() != s2.size())
        return false;
    int res = 0;
    for (int i = 0; i < s1.size(); ++i)
    {
        if (s1[i] != s2[i])
            ++res;
        if (res > d)
            return false;
    }
    return true;
}

class kmerProcessor
{
    int _k;
    string pat;
    string const alphabet = "ACGT";
    std::ofstream ost;
    void processAllKmersRec(string const&, int const, int const);
public:
    kmerProcessor(string const& pattern) : _k(pattern.size()), pat(pattern)
    {
        ost = std::ofstream(OUTFILENAME);
        if (!ost) 
        {
        std::cerr << "Cannot open output file!\n";
        exit(1);
        }
    }
    void processAllKmers(int const);
};

void kmerProcessor::processAllKmersRec(string const& prefix, int const k, int const d)
{
    if (k == 0) // resursion stop: process current kmer
    {   
        if (hammingDistanceLEd(prefix, pat, d)) // if dist <= d
            ost << prefix << '\n';
        return;
    }
    for (int i = 0; i < alphabet.size(); ++i)
    {   // far all possible combinations in ith pos
        string newPrefix;
        newPrefix = prefix + alphabet[i];
        processAllKmersRec(newPrefix, k-1, d);
    }
}

void kmerProcessor::processAllKmers(int const d)
{
    processAllKmersRec("", _k, d);
}

int main()
{   
    std::ifstream ist{INFILENAME};
    if (!ist) 
    {
        std::cerr << "Cannot open input file!\n";
        exit(1);
    }

    string text;
    int d;
    getline(ist, text);
    ist >>  d;

    kmerProcessor processor(text);
    processor.processAllKmers(d);
    return 0;
}