#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>

using std::vector;
using std::string;
using std::set;
const char* INFILENAME = "rosalind_ba2b.txt";
const char* OUTFILENAME = "out.txt";

int hammingDistance(string const& s1, string const& s2)
{
    int res = 0;
    for (int i = 0; i < s1.size(); ++i)
        if (s1[i] != s2[i])
            ++res;
    return res;
}

class kmerProcessor
{
    int K;
    int min_d;
    string result;
    string const alphabet = "ACGT";
    vector<set<string>> dnaPatterns;
    void addPatterns(string const&);
    int minDistance(string const&, int const);
    void processAllKmersRec(string const&, int const);
public:
    kmerProcessor(std::ifstream&);
    void processAllKmers(std::ofstream&);
};

kmerProcessor::kmerProcessor(std::ifstream& ist)
{
    ist >> K;
    ist.ignore(1, '\n');
    string tmp;
    while (getline(ist, tmp))
        addPatterns(tmp);
    min_d = INT32_MAX;
}

void kmerProcessor::addPatterns(string const& text)
{
    set<string> patterns;
    for (int i = 0; i <= (text.size()-K); ++i)
        patterns.insert(text.substr(i, K));
    dnaPatterns.push_back(patterns);
}

int kmerProcessor::minDistance(string const& kmer, int const idx)
{
    int result = INT32_MAX;
    for (auto pat: dnaPatterns[idx])
    {
        int dist = hammingDistance(kmer, pat);
        if (dist < result)
            result = dist;
    }
    return result;
}

void kmerProcessor::processAllKmersRec(string const& prefix, int const k)
{
    if (k == 0) // resursion stop: process current kmer
    {   
        int cur_dist = 0;
        for (int i = 0; i < dnaPatterns.size(); ++i)
            cur_dist += minDistance(prefix, i);
        if (cur_dist < min_d)
        {
            min_d = cur_dist;
            result = prefix;
        }
        return;
    }
    for (int i = 0; i < alphabet.size(); ++i)
    {   // far all possible combinations in ith pos
        string newPrefix;
        newPrefix = prefix + alphabet[i];
        processAllKmersRec(newPrefix, k-1);
    }
}

void kmerProcessor::processAllKmers(std::ofstream& ost)
{
    processAllKmersRec("", K);
    ost << result << '\n';
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
    processor.processAllKmers(ost);
    return 0;
}