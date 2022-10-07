#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>

using std::set;
using std::string;
using std::vector;
const char* INFILENAME = "rosalind_ba2a.txt";
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
    int k, d;
    set<string> kmers;
    set<string> patterns;
    vector<string> dna;
    set<string> genNeighbors(string const&);
    bool containsApproxMatch(string const&, string const&);
public:
    kmerProcessor(std::ifstream& ist)
    {
        ist >> k >> d;
        ist.ignore(INT32_MAX, '\n');
        string tmp;
        while (getline(ist, tmp))
        {
            for (int i = 0; i <= (tmp.size()-k); ++i)
                kmers.insert(tmp.substr(i, k));
            dna.push_back(tmp);
        }
    }
    void MotifEnumeration();
    set<string> getPatterns() const
    {
        return patterns;
    }
};
// recursive function to get all less-than-d neighbors
set<string> kmerProcessor::genNeighbors(string const& kmer)
{
    set<string> neighbors;
    static const vector<string> BASES = {"A", "C", "G", "T"};
    if (d == 0)
    {
        neighbors.insert(kmer);
        return neighbors;
    }
    if (kmer.size() == 1)
    {
        for (string item: BASES)
            neighbors.insert(item);
        return neighbors;
    }
    for (auto nei: genNeighbors(kmer.substr(1)))
    {
        if (hammingDistanceLEd(kmer.substr(1), nei, d-1))
        {
            for (string base: BASES)
                neighbors.insert(base + nei);
        } else {
            neighbors.insert(kmer.substr(0, 1) + nei);
        }
    }
    return neighbors;
}

bool kmerProcessor::containsApproxMatch(string const& s, string const& kmer)
{

    for (int i = 0; i <= (s.size()-k); ++i)
    {
        if (hammingDistanceLEd(s.substr(i, k), kmer, d))
            return true;
    }
    return false;
}

void kmerProcessor::MotifEnumeration()
{
    patterns.clear();
    for (string kmer: kmers)
    {
        for (string nei: genNeighbors(kmer))
        {
            bool isMatched = true;
            for (string seq: dna)
            {
                isMatched &= containsApproxMatch(seq, nei);
                if (!isMatched)
                    break;
            }
            if (isMatched)
                patterns.insert(nei);
        }
    }
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
    processor.MotifEnumeration();
    for (auto pat: processor.getPatterns())
        ost << pat << " ";
    ost << '\n';
    return 0;
}