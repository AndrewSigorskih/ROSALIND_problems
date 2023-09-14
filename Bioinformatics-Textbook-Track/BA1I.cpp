#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

using std::map;
using std::string;
using std::vector;
const char* INFILENAME = "rosalind_ba1i.txt";
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
    int max_freq; // maximum observed frequency
    string const alphabet = "ACGT";
    map<string, int> kmers;
    map<int, vector<string>> dict; // stores groups of kmers
    void processAllKmersRec(string const&, int const, int const);
public:
    kmerProcessor(string const& text, int k)
    {
        _k = k;
        max_freq = 0;
        for (int i = 0; i <= (text.size()-k); ++i)
            kmers[text.substr(i, k)]++;
    }
    void printKmers(std::ofstream&);
    void processAllKmers(int const, std::ofstream&);
};

void kmerProcessor::printKmers(std::ofstream& ost)
{
    for (auto const& item: kmers)
        ost << item.first << " ";
    ost << "\n";
}

void kmerProcessor::processAllKmersRec(string const& prefix, int const k, int const d)
{
    if (k == 0) // resursion stop: process current kmer
    {   
        int count = 0; // how many kmers from text match current one
        for (auto const& item: kmers)
        {
            if (hammingDistanceLEd(prefix, item.first, d)) // if dist <= d
            {
                count += item.second;
            }
        }
        if (count > 0) // store only good-scoring kmers
            dict[count].push_back(prefix);
        if (count > max_freq) // store best count
            max_freq = count;
        return;
    }
    for (int i = 0; i < alphabet.size(); ++i)
    {   // far all possible combinations in ith pos
        string newPrefix;
        newPrefix = prefix + alphabet[i];
        processAllKmersRec(newPrefix, k-1, d);
    }
}

void kmerProcessor::processAllKmers(int const d, std::ofstream& ost)
{
    processAllKmersRec("", _k, d);
    for (auto item: dict[max_freq])
        ost << item << " ";
    ost << '\n';
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

    string text;
    int k, d;
    getline(ist, text);
    ist >> k >> d;

    kmerProcessor processor(text, k);
    processor.processAllKmers(d, ost);
    return 0;
}